/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/integral.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/lib3index/dfhelper.h"

#include "jk.h"

#include <sstream>
#include "psi4/libpsi4util/PsiOutStream.h"
#ifdef _OPENMP
#include <omp.h>
#include "psi4/libpsi4util/process.h"
#endif

using namespace psi;

namespace psi {
/* compute_D needs to be run before this functions. */
/*Mem_2B_DFJK::allocate_JK() {
    compute_D();
}*/
Mem_2B_DFJK::Mem_2B_DFJK( std::shared_ptr<BasisSet> primary, 
                        std::shared_ptr<BasisSet> auxiliary) 
                        : JK(primary), 
                          auxiliary_(auxiliary) {
    common_init();
}
Mem_2B_DFJK::Mem_2B_DFJK( std::shared_ptr<BasisSet> primary, 
                        std::shared_ptr<BasisSet> auxiliary,
                        std::shared_ptr<BasisSet> one_bas,
                        std::shared_ptr<BasisSet> two_bas) 
                        : JK(primary), 
                          auxiliary_(auxiliary),
                          one_basis_(one_bas), 
                          two_basis_(two_bas)
                          {
    common_init(true);
    nbf_ = primary_->nbf();
    naux_ = auxiliary_->nbf();
    nob_ = one_basis_->nbf();
    ntb_ = two_basis_->nbf();
}
Mem_2B_DFJK::~Mem_2B_DFJK() {}
void Mem_2B_DFJK::common_init() { 
    outfile->Printf("Hello from common_init()\n");
    dfh_ = std::make_shared<DFHelper>(primary_, auxiliary_); 
}
void Mem_2B_DFJK::common_init(bool k) { dfh_ = std::make_shared<DFHelper_2B>(one_basis_, two_basis_, primary_, auxiliary_); }
size_t Mem_2B_DFJK::memory_estimate() {
    dfh_->set_nthreads( df_ints_num_threads_ );
    dfh_->set_schwarz_cutoff(cutoff_);
//joejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoe
// we need to change this function for when secondary is set
    return dfh_->get_core_size();
//joejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoe
}
void Mem_2B_DFJK::preiterations() {
    // Pass our information to DFHelper
    dfh_->set_nthreads(omp_nthread_);
    dfh_->set_schwarz_cutoff(cutoff_);
    dfh_->set_method("STORE");
    dfh_->set_fitting_condition(condition_);
    dfh_->set_memory(memory_ - memory_overhead());
    dfh_->set_do_wK(do_wK_);
    dfh_->set_omega(omega_);

    // DFHelper prepares the quantities it needs for SCF iterations
    dfh_->initialize();
}
void Mem_2B_DFJK::compute_D_2B(){

    if (C_left_.size() && !C_right_.size() ) {
        lr_symmetric_ = true;
    } else {
        lr_symmetric_ = false;
    }

    /* as I started writing this function I realized I needed more functionality
       than the name suggests. This function does as follows:
       Prepares shared matrices for coulomb and exchange partitions
       Prepares SCF Density matrices
       Copies coefficient and density matrices into their AO versions.*/
    J_oo_.clear(); /*You may wonder, why clear these if we aren't goint to use */
    J_ot_.clear(); /* them right now? We want to actively make sure that the ao*/
    J_tt_.clear(); /* coulomb and exchange matrices get used*/
    J_oo_ao_.clear();
    J_ot_ao_.clear();
    J_tt_ao_.clear();
    
    K_oo_.clear();
    K_ot_.clear();
    K_tt_.clear();
    K_oo_ao_.clear();
    K_ot_ao_.clear();
    K_tt_ao_.clear();

    D_ao_.clear();

    //size_t nocc = static_cast<size_t>(max_nocc());
    C_left_ao_ = C_left_;
    if (lr_symmetric_) {
        C_right_ = C_left_;
        C_right_ao_ = C_left_ao_;
    } else { /* it's on the user to provide C_left in this case */
        C_right_ao_ = C_right_;
    }
    if (C_left_ao_.size() != C_right_ao_.size()) {
        std::stringstream error;
        error << "Mem_2B_DFJK::compute_D_2B():C_right_.size() doesn't match \nC_left_.size()!\n";
        throw PSIEXCEPTION(error.str().c_str());
    }

    size_t maxx_nocc = max_nocc();

    /* I know I called this function compute_D, but computing D comes with making
       sure we have containers for all the data we actually plan to store. */
    /* D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_*/
    for (int i = 0; i < C_left_ao_.size(); i++){
        double* clp = C_left_ao_[i]->pointer()[0];
        double* crp = C_right_ao_[i]->pointer()[0];
        int cl_nocc = C_left_ao_[i]->coldim(0);
        int cr_nocc = C_right_ao_[i]->coldim(0);
        std::stringstream Ds;
        Ds << "D " << i << " (AO)";
        D_ao_.push_back(std::make_shared<Matrix>(Ds.str(), nbf_, nbf_));
        double* dp = D_ao_[i]->pointer()[0];
        C_DGEMM('N', 'T', nbf_, nbf_, maxx_nocc, 1.0, clp, cl_nocc, crp, cr_nocc, 0.0, dp, nbf_);
    }
    /* D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_D_ao_*/
    for (int i = 0; i < C_left_ao_.size(); i++) {
        std::stringstream s;
        s << "J oo " << i << " (AO)";
        J_oo_ao_.push_back(std::make_shared<Matrix>(s.str(), nob_, nob_));
    }
    for (int i = 0; i < C_left_ao_.size(); i++) {
        /* for our initial application, I could just do nbf instead of ntb, but
           I don't want to make anyone have to come back and fix it later, so
           we're going to use ntb_ in this function. 
           Let the user be warned, if you're including the AO basis, you should
           be including it as primary_ and as two_basis_*/
        std::stringstream s;
        s << "J ot " << i << " (AO)";
        J_ot_ao_.push_back(std::make_shared<Matrix>(s.str(), nob_, ntb_));
    }
    for (int i = 0; i < C_left_ao_.size(); i++) {
        std::stringstream s;
        s << "J tt " << i << " (AO)";
        J_tt_ao_.push_back(std::make_shared<Matrix>(s.str(), ntb_, ntb_));
    }
    for (int i = 0; i < C_left_ao_.size(); i++) {
        std::stringstream s;
        s << "K oo " << i << " (AO)";
        K_oo_ao_.push_back(std::make_shared<Matrix>(s.str(), nob_, nob_));
    }
    for (int i = 0; i < C_left_ao_.size(); i++) {
        std::stringstream s;
        s << "K ot " << i << " (AO)";
        K_ot_ao_.push_back(std::make_shared<Matrix>(s.str(), nob_, ntb_));
    }
    for (int i = 0; i < C_left_ao_.size(); i++) {
        std::stringstream s;
        s << "K tt " << i << " (AO)";
        K_tt_ao_.push_back(std::make_shared<Matrix>(s.str(), ntb_, ntb_));
    }
}
    /* This functions exists because it is necessary for Psi4 to compile
     *     works the same as MemDFJK::compute_JK()
     */
void Mem_2B_DFJK::compute_JK() {
    dfh_->build_JK(C_left_ao_, C_right_ao_, D_ao_, J_ao_, K_ao_, wK_ao_, max_nocc(), do_J_, do_K_, do_wK_, lr_symmetric_);
    if (lr_symmetric_) {
        if (do_wK_) {
            for (size_t N = 0; N < wK_ao_.size(); N++) {
                wK_ao_[N]->hermitivitize();
            }
        }
    }
}
/* We're going to assume that twob is the AO basis set*/
void Mem_2B_DFJK::compute_2B_JK() {
    /* for now, we're going to boldly assume (require that our user ensure) 
       that our system is C1*/
    compute_D_2B();
    dfh_->build_2B_JK(C_left_ao_, C_right_ao_, D_ao_, J_tt_ao_, J_ot_ao_, J_oo_ao_, K_tt_ao_, K_ot_ao_, K_oo_ao_, static_cast<size_t>(max_nocc()), do_J_, do_K_, lr_symmetric_);
}
void Mem_2B_DFJK::postiterations() {}
void Mem_2B_DFJK::print_header() const {
    if (print_) {
        outfile->Printf("  ==> Mem_2B_DFJK: Density-Fitted J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:           %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    K tasked:           %11s\n", (do_K_ ? "Yes" : "No"));
        outfile->Printf("    wK tasked:          %11s\n", (do_wK_ ? "Yes" : "No"));
        if (do_wK_) outfile->Printf("    Omega:              %11.3E\n", omega_);
        outfile->Printf("    OpenMP threads:     %11d\n", omp_nthread_);
        outfile->Printf("    Memory [MiB]:       %11ld\n", (memory_ * 8L) / (1024L * 1024L));
        outfile->Printf("    Algorithm:          %11s\n", (dfh_->get_AO_core() ? "Core" : "Disk"));
        outfile->Printf("    Schwarz Cutoff:     %11.0E\n", cutoff_);
        outfile->Printf("    Mask sparsity (%%):  %11.4f\n", 100. * dfh_->ao_sparsity());
        outfile->Printf("    Fitting Condition:  %11.0E\n\n", condition_);

        outfile->Printf("   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
    }
}
int Mem_2B_DFJK::max_nocc() const {
    int max_noc = 0;
    for (size_t N = 0; N < C_left_ao_.size(); N++) {
        max_noc = ( C_left_ao_[N]->coldim(0) > max_noc ? C_left_ao_[N]->coldim(0) : max_noc); 
    }
    return max_noc;
}
void Mem_2B_DFJK::set_do_wK(bool tf) { do_wK_ = tf; dfh_->set_do_wK(tf); }
void Mem_2B_DFJK::set_do_JK_oo(bool v) {do_JK_oo_=v; if (dfh_) {dfh_->set_do_JK_oo(v);}}
void Mem_2B_DFJK::set_do_JK_ot(bool v) {do_JK_ot_=v; if (dfh_) {dfh_->set_do_JK_ot(v);}}
void Mem_2B_DFJK::set_do_JK_tt(bool v) {do_JK_tt_=v; if (dfh_) {dfh_->set_do_JK_tt(v);}}
} // namespace psi
