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

MemDF_2B_JK::MemDF_2B_JK( std::shared_ptr<BasisSet> primary, 
                        std::shared_ptr<BasisSet> auxiliary) 
                        : JK(primary), 
                          auxiliary_(auxiliary) {
    common_init();
}

MemDF_2B_JK::MemDF_2B_JK( std::shared_ptr<BasisSet> primary, 
                        std::shared_ptr<BasisSet> row_bas,
                        std::shared_ptr<BasisSet> col_bas, 
                        std::shared_ptr<BasisSet> auxiliary) 
                        : JK(primary), 
                          row_bas_(row_bas), 
                          col_bas_(col_bas),
                          auxiliary_(auxiliary) {
    common_init(true);
    printf("ping\n");
}
MemDF_2B_JK::~MemDF_2B_JK() {}
void MemDF_2B_JK::common_init() { dfh_ = std::make_shared<DFHelper>(primary_, auxiliary_); }
void MemDF_2B_JK::common_init(bool k) { dfh_ = std::make_shared<DFHelper>(primary_, row_bas_, col_bas_, auxiliary_); }
size_t MemDF_2B_JK::memory_estimate() {
    dfh_->set_nthreads( df_ints_num_threads_ );
    dfh_->set_schwarz_cutoff(cutoff_);
//joejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoe
// we need to change this function for when secondary is set
    return dfh_->get_core_size();
//joejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoejoe
}

void MemDF_2B_JK::preiterations() {
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
    /* This functions exists because it is necessary for Psi4 to compile
     *     works the same as MemDFJK::compute_JK()
     */
void MemDF_2B_JK::compute_JK() {
    dfh_->build_JK(C_left_ao_, C_right_ao_, D_ao_, J_ao_, K_ao_, wK_ao_, max_nocc(), do_J_, do_K_, do_wK_, lr_symmetric_);
    if (lr_symmetric_) {
        if (do_wK_) {
            for (size_t N = 0; N < wK_ao_.size(); N++) {
                wK_ao_[N]->hermitivitize();
            }
        }
    }

}
void MemDF_2B_JK::postiterations() {}
void MemDF_2B_JK::print_header() const {
    if (print_) {
        outfile->Printf("  ==> MemDFJK: Density-Fitted J/K Matrices <==\n\n");

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
int MemDF_2B_JK::max_nocc() const {
    int max_nocc = 0;
    for (size_t N = 0; N < C_left_ao_.size(); N++) {
        max_nocc = (C_left_ao_[N]->colspi()[0] > max_nocc ? C_left_ao_[N]->colspi()[0] : max_nocc); 
    }
    max_nocc;
}
void MemDF_2B_JK::set_do_wK(bool tf) { do_wK_ = tf; dfh_->set_do_wK(tf); }
} // namespace psi

