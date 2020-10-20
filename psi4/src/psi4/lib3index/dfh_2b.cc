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

#include "dfhelper.h"

#include <algorithm>
#include <cstdlib>
#ifdef _MSC_VER
#include <process.h>
#define SYSTEM_GETPID ::_getpid
#else
#include <unistd.h>
#define SYSTEM_GETPID ::getpid
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/psi4-dec.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libfock/jk.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/sieve.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"

#include "dftensor.h"

namespace psi {

/* The Basis Sets should suffice to specify the object */
DFHelper_2B::DFHelper_2B( std::shared_ptr<BasisSet> one_basis, std::shared_ptr<BasisSet> two_basis, std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> aux ) : DFHelper(primary, aux) {
    one_basis_ = one_basis;
    two_basis_ = two_basis;

    nob_ = one_basis_->nbf();
    ntb_ = two_basis_->nbf();

    oshells_ = one_basis_->nshell();
    tshells_ = two_basis_->nshell();
    prepare_blocking();
    form_prim_to_ob();
}

//DFHelper_2B::~DFHelper() {}

// useful for f12 sapt. We need to be able to get the same functions from
//   the primary basis in the cabs basis

void DFHelper_2B::form_prim_to_ob() {
    prim_to_ob_.resize(nbf_);
    prim_to_ob_sh_.resize(pshells_);
    // number of atoms in the basis set
    std::vector<size_t> ob_atom_to_shell(primary_->shell(primary_->nshell()-1).ncenter()+1);
    int lastatom = -1;
    for (int i = 0; i < one_basis_->nshell(); i++) {
        if (one_basis_->shell(i).ncenter() != lastatom ) {
            lastatom = one_basis_->shell(i).ncenter();
            ob_atom_to_shell[lastatom] = i;
        }
    }
    int shell_on_atom_index = -1;
    int atom_last = -1;
    for (int i = 0; i < primary_->nshell(); i++) {
        if (primary_->shell(i).ncenter() != atom_last) {
            atom_last = primary_->shell(i).ncenter();
            shell_on_atom_index = 0;
        } else {
            shell_on_atom_index++;
        }
        int pr_shell_start = primary_->shell(i).function_index();
        int ob_shell_start = one_basis_->shell(
ob_atom_to_shell[atom_last] +shell_on_atom_index
                                               ).function_index();
        prim_to_ob_sh_[i] = ob_atom_to_shell[atom_last] + shell_on_atom_index;
        for (int j = 0; j < primary_->shell(i).nfunction(); j++ ) {
            prim_to_ob_[pr_shell_start + j] = ob_shell_start + j;
        }
    }
    ob_to_prim_.resize(nob_, -1ul);
    ob_to_prim_sh_.resize(oshells_, -1ul);
    for (size_t i = 0; i < nbf_; i++) {
        ob_to_prim_[ prim_to_ob_[i] ] = i;
    }
    for (size_t i = 0; i < pshells_; i++) {
        ob_to_prim_sh_[ prim_to_ob_sh_[i] ] = i;
    }
    lastind_.resize(nob_);
    size_t lastnum;
    for (size_t i = 0; i < nob_; i++) {
        lastind_[i] = (ob_to_prim_[i] == -1 ? lastnum : ob_to_prim_[i] );
        lastnum = (ob_to_prim_[i] == -1 ? lastnum : ob_to_prim_[i] );
    }
}

void DFHelper_2B::initialize() {

    hold_met_ = true;
    /* This is a provisional step */
    AO_core();


    if (!(std::fabs(mpower_ - 0.0) < 1e-13)) (hold_met_ ? prepare_metric_core() : prepare_metric());
    prepare_metric_lu_core();

    prepare_sparsity();

    if (two_basis_.get() == primary_.get()) {
        if (AO_core_) {
            prepare_AO_core();
        } else {
            std::stringstream error;
            error << "DFHelper_2B:initialize:Out of core machinery not written\n";
            throw PSIEXCEPTION(error.str().c_str());
        }
    } else {
        std::stringstream error;
        error << "DFHelper_2B:initialize: one_basis must be the same as AO_basis. Generic Case Currently Unsupported\n";
        throw PSIEXCEPTION(error.str().c_str());
    }
}

void DFHelper_2B::AO_core() {
    prepare_sparsity();
    required_core_size_  = big_skips_[nbf_];
    required_core_size_ += ob_ao_big_skips_[nob_];
    required_core_size_ += ob_ob_big_skips_[nob_];
    required_core_size_ += 2 * naux_ * naux_;
    required_core_size_ += nthreads_ * nbf_ * nbf_;
    required_core_size_ += 3 * nob_ * nob_ * Qshell_max_;
    AO_core_ = required_core_size_ < memory_;
}

size_t DFHelper_2B::calcfull_3index(bool symm) {
    return big_skips_[nbf_] + ob_ao_big_skips_[nob_] + ob_ob_big_skips_[nob_];
}

void DFHelper_2B::prepare_blocking() {
    Qshells_ = aux_->nshell();
    pshells_ = primary_->nshell();
    oshells_ = one_basis_->nshell();
    tshells_ = two_basis_->nshell();

    Qshell_aggs_.resize(Qshells_ + 1);
    pshell_aggs_.resize(pshells_ + 1);
    oshell_aggs_.resize(oshells_ + 1);
    tshell_aggs_.resize(tshells_ + 1);

    Qshell_max_ = aux_->max_function_per_shell();
    // AUX shell blocking 
    Qshell_aggs_[0] = 0;
    for (size_t i = 0; i < Qshells_; i++) {
        Qshell_aggs_[i + 1] = Qshell_aggs_[i] + aux_->shell(i).nfunction();
    }

    // AO shell blocking
    pshell_aggs_[0] = 0;
    for (size_t i = 0; i < pshells_; i++) {
        pshell_aggs_[i + 1] = pshell_aggs_[i] + primary_->shell(i).nfunction();
    }

    // AO shell blocking
    oshell_aggs_[0] = 0;
    for (size_t i = 0; i < oshells_; i++) {
        oshell_aggs_[i + 1] = oshell_aggs_[i] + one_basis_->shell(i).nfunction();
    }

    // AO shell blocking
    tshell_aggs_[0] = 0;
    for (size_t i = 0; i < tshells_; i++) {
        tshell_aggs_[i + 1] = tshell_aggs_[i] + two_basis_->shell(i).nfunction();
    }
}

void DFHelper_2B::prepare_sparsity() {
    if (sparsity_prepared_) {return;}

    timer_on("DFH_2B: sparsity prep");
    double max_integral = 0.0;
    double ob_ao_max_integral = 0.0;
    double ob_ob_max_integral = 0.0;

    std::vector<double> shell_max_vals(pshells_ * pshells_ , 0.0);
    std::vector<double> ob_ao_shell_max_vals(pshells_ * oshells_, 0.0);
    std::vector<double> ob_ob_shell_max_vals(oshells_ * oshells_, 0.0);
    std::vector<double> fun_max_vals(nbf_ * nbf_, 0.0);
    std::vector<double> ob_ao_fun_max_vals(nbf_ * nob_, 0.0);
    std::vector<double> ob_ob_fun_max_vals(nob_ * nob_, 0.0);

    small_skips_.resize(nbf_ + 1);
    big_skips_.resize(nbf_ + 1);
    symm_small_skips_.resize(nbf_);
    symm_big_skips_.resize(nbf_ + 1);
    symm_ignored_columns_.resize(nbf_);


    ob_ao_small_skips_.resize(nob_ + 1);
    ob_ao_big_skips_.resize(nob_ + 1);

    ob_ob_small_skips_.resize(nob_ + 1);
    ob_ob_big_skips_.resize(nob_ + 1);

    ob_ao_small_skips_[0] = nbf_;
    ob_ao_big_skips_[0] = 0;

    ob_ob_small_skips_[0] = nob_;
    ob_ob_big_skips_[0] = 0;

    schwarz_fun_mask_.resize(nbf_ * nbf_);
    schwarz_shell_mask_.resize(pshells_ * pshells_);

    ob_ao_schwarz_fun_mask_.resize(nob_ * nbf_);//oshells_ * pshells_);
    ob_ao_schwarz_shell_mask_.resize(oshells_ * pshells_);

    ob_ob_schwarz_fun_mask_.resize(nob_ * nob_);
    ob_ob_schwarz_shell_mask_.resize(nob_ * nob_);

    size_t nthreads = (nthreads_ == 1 ? 1 : 2);
    auto rifactory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);
    auto ob_ao_rifactory = std::make_shared<IntegralFactory>(one_basis_, primary_, one_basis_, primary_);
    auto ob_ob_rifactory = std::make_shared<IntegralFactory>(one_basis_, one_basis_, one_basis_, one_basis_);

    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(nthreads);
    std::vector<std::shared_ptr<TwoBodyAOInt>> ob_ao_eri(nthreads);
    std::vector<std::shared_ptr<TwoBodyAOInt>> ob_ob_eri(nthreads);
    std::vector<const double*> buffer(nthreads);
    std::vector<const double*> ob_ao_buffer(nthreads);
    std::vector<const double*> ob_ob_buffer(nthreads);


#pragma omp parallel num_threads(nthreads) if (nbf_ > 1000)
{
    int rank = 0;
#ifdef _OPENMP
    rank = omp_get_thread_num();
#endif
    eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    ob_ao_eri[rank] = std::shared_ptr<TwoBodyAOInt>(ob_ao_rifactory->eri());
    ob_ob_eri[rank] = std::shared_ptr<TwoBodyAOInt>(ob_ob_rifactory->eri());

    buffer[rank] = eri[rank]->buffer();
    ob_ao_buffer[rank] = ob_ao_eri[rank]->buffer();
    ob_ob_buffer[rank] = ob_ob_eri[rank]->buffer();
}


//#pragma omp parallel for num_threads(nthreads)
    for (size_t ob_iters = 0; ob_iters < oshells_; ob_iters++) {
        int rank = 0;
//#ifdef _OPENMP
        rank = omp_get_thread_num();
//#endif
        size_t nphi = one_basis_->shell(ob_iters).nfunction();
        size_t ind_phi = one_basis_->shell(ob_iters).function_index();
        for (size_t ob_jters = 0; ob_jters < oshells_; ob_jters++) {
            ob_ob_eri[rank]->compute_shell(ob_iters, ob_jters, ob_iters, ob_jters);
            size_t ngam = one_basis_->shell(ob_jters).nfunction();
            size_t ind_gamma = one_basis_->shell(ob_jters).function_index();
            for (size_t ob_iter = 0; ob_iter < nphi; ob_iter++) {
                size_t phi = ind_phi + ob_iter;
                for (size_t ob_jter = 0; ob_jter < ngam; ob_jter++) {
                    size_t gamma = ind_gamma + ob_jter;
                    double val = fabs(ob_ob_buffer[rank][ ob_iter * (ngam + ngam * nphi * ngam) + ob_jter * (1 + ngam * nphi ) ]);
                    ob_ob_max_integral = std::max(ob_ob_max_integral, val);
                    ob_ob_shell_max_vals[ ob_iters * oshells_ + ob_jters ] = std::max(ob_ob_shell_max_vals[ ob_iters * oshells_ + ob_jters ], val);
                    ob_ob_fun_max_vals[ phi * nob_ + gamma ] = std::max(ob_ob_fun_max_vals[ phi * nob_ + gamma ], val);
                    if ( ob_to_prim_[gamma] != -1ul ) {
                        ob_ao_fun_max_vals[ phi * nbf_ + ob_to_prim_[gamma] ] = std::max(ob_ao_fun_max_vals[ phi * nbf_ + ob_to_prim_[gamma] ], val);
                        if (ob_to_prim_[phi] != -1ul) {
                            fun_max_vals[ ob_to_prim_[phi] * nbf_ + ob_to_prim_[gamma]] = std::max(fun_max_vals[ ob_to_prim_[phi] * nbf_ + ob_to_prim_[gamma]], val);
                        }
                    }
                    if ( ob_to_prim_sh_[ob_jters] != -1ul ) {
                        ob_ao_shell_max_vals[ob_iters * pshells_ + ob_to_prim_sh_[ob_jters] ] = std::max(ob_ao_shell_max_vals[ob_iters * pshells_ + ob_to_prim_sh_[ob_jters] ], val);
                        if (ob_to_prim_sh_[ob_iters] != -1ul){
                            shell_max_vals[ ob_to_prim_sh_[ob_iters] * pshells_ + ob_to_prim_sh_[ob_jters]] = std::max(shell_max_vals[ ob_to_prim_sh_[ob_iters] * pshells_ + ob_to_prim_sh_[ob_jters]], val);
                        }
                    }
                }
            }
        }
    }

    double tolerance = cutoff_ * cutoff_/ ob_ob_max_integral;

    for (size_t ao_iters = 0; ao_iters < pshells_ * pshells_; ao_iters++) {
        schwarz_shell_mask_[ao_iters] = (shell_max_vals[ao_iters] < tolerance ? 0 : 1);
    }

    for (size_t ao_iter = 0, count = 0; ao_iter < nbf_; ao_iter++) {
        count = 0;
        for (size_t ao_jter = 0; ao_jter < nbf_; ao_jter++) {
            if (fun_max_vals[ao_iter * nbf_ + ao_jter] > tolerance ) {
                count++;
                schwarz_fun_mask_[ao_iter * nbf_ + ao_jter] = count;
            } else {
                schwarz_fun_mask_[ao_iter * nbf_ + ao_jter] = 0;
            }
        }
        small_skips_[ao_iter] = count;
    }

    big_skips_[0] = 0;
    size_t coltots = 0; 
    for (size_t ao_iter = 1; ao_iter < nbf_ + 1; ao_iter++) {
        coltots += small_skips_[ao_iter - 1];
        big_skips_[ao_iter] = big_skips_[ao_iter - 1] + naux_ * small_skips_[ao_iter - 1];
    }
    small_skips_[nbf_] = coltots;

    // thanks matt
    for (size_t i = 0; i < nbf_; i++) {
        size_t size = 0;
        size_t skip = 0;
        for (size_t j = 0; j < nbf_; j++) {
            if (schwarz_fun_mask_[i * nbf_ + j]) {
                (j >= i ? size++: skip++);
            }
        }
        symm_small_skips_[i] = size;
        symm_ignored_columns_[i] = skip;
    }

    for (size_t ob_iters = 0; ob_iters < oshells_ * pshells_; ob_iters++ ) {ob_ao_schwarz_shell_mask_[ob_iters] = ( ob_ao_shell_max_vals[ob_iters] > tolerance ? 1 : 0); } 

    for (size_t ob_iter = 0, count = 0; ob_iter < nob_; ob_iter++) {
        count = 0;
        for (size_t ao_iter = 0; ao_iter < nbf_; ao_iter++) {
            if (ob_ao_fun_max_vals[ob_iter * nbf_ + ao_iter] > tolerance) {
                count++;
                ob_ao_schwarz_fun_mask_[ob_iter * nbf_ + ao_iter] = count;
            } else {
                ob_ao_schwarz_fun_mask_[ob_iter * nbf_ + ao_iter] = 0;
            }
        }
        ob_ao_small_skips_[ob_iter] = count;
    } 

    coltots = 0;
    ob_ao_big_skips_[0] = 0;
    for (size_t ob_iter = 0; ob_iter < nob_; ob_iter++) {
        ob_ao_big_skips_[ob_iter + 1] = ob_ao_big_skips_[ob_iter] + ob_ao_small_skips_[ob_iter] * naux_;
        coltots += ob_ao_small_skips_[ob_iter];
    }
    ob_ao_small_skips_[nob_] = coltots;


    for (size_t ob_iters = 0; ob_iters < oshells_ * oshells_; ob_iters++) { ob_ob_schwarz_shell_mask_[ob_iters] = ( ob_ob_shell_max_vals[ob_iters] > tolerance ? 1 : 0); }


    for (size_t ob_iter = 0, count = 0; ob_iter < nob_; ob_iter++) {
        count = 0;
        for (size_t ob_jter = 0; ob_jter < nob_; ob_jter++) {
            if (ob_ob_fun_max_vals[ob_iter * nob_ + ob_jter] > tolerance) {
                count++;
                ob_ob_schwarz_fun_mask_[ ob_iter * nob_ + ob_jter ] = count;
            } else {
                ob_ob_schwarz_fun_mask_[ ob_iter * nob_ + ob_jter ] = 0;
            }
        }
        ob_ob_small_skips_[ob_iter] = count;
    }

    ob_ob_big_skips_[0] = 0;
    coltots = 0;
    for (size_t ob_iter = 0; ob_iter < nob_; ob_iter++ ) {
        ob_ob_big_skips_[ob_iter+ 1] = ob_ob_big_skips_[ob_iter] + ob_ob_small_skips_[ob_iter] * naux_;
        coltots += ob_ob_small_skips_[ob_iter] ;
    }
    ob_ob_small_skips_[nob_] = coltots;

    sparsity_prepared_ = true;
    timer_off("DFH_2B: sparsity prep");
} // DFHelper_2B::prepare_sparsity() 

std::pair<size_t, size_t> DFHelper_2B::fshell_blocks_for_ob_ao_AO_build(const size_t mem,
                                                     size_t symm,
                                    std::vector<std::pair<size_t, size_t>>& b) {
    size_t full_3index = calcfull_3index(symm);
    size_t demand, end, begin, current, block_size, tmpbs, total, count, largest;
    block_size = tmpbs= total= count= largest = 0;
    for (size_t i = 0; i < oshells_; i++) {
        count++;
        begin = oshell_aggs_[i];
        end = oshell_aggs_[i + 1] - 1;
        tmpbs += end - begin + 1; // I think temporary basis functions
        if (symm) {
           std::stringstream error;
           error << "Don't use Symm code for two-basis AOs.";
           throw PSIEXCEPTION(error.str().c_str()); 
        } else {
            current = ob_ao_big_skips_[end+1] - ob_ao_big_skips_[begin];
            total += 2*current;
        }

        demand = total;
        demand += full_3index;
        demand += (hold_met_ ? naux_ * naux_ : total);
        if ( demand > mem || i == oshells_ -1 ) {
            if ( count == 1 && i != oshells_ - 1 ) {
                std::stringstream error;
                error << "DFHelper: not enough memory for (p shell) AO blocking!"
                      << " required memory: " << demand * 8 / (1024 * 1024 * 1024.0) << " [GiB].";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if ( demand > mem ) {
                total -= current;
                tmpbs -= end - begin + 1;
                b.push_back(std::make_pair(i - count + 1, i - 1));
                i--;
            } else if (i == oshells_ - 1) {
                b.push_back(std::make_pair(i-count+1, i));
            }
            if (largest < total) {
                largest = total;
                block_size = tmpbs;
            }
            count = 0;
            total = 0;
            tmpbs = 0;
        }
    }
    return std::make_pair(largest, block_size);
}
std::pair<size_t, size_t> DFHelper_2B::fshell_blocks_for_ob_ob_AO_build(const size_t mem,
                                                     size_t symm,
                                    std::vector<std::pair<size_t, size_t>>& b) {
    size_t full_3index = calcfull_3index(symm);
    size_t demand, end, begin, current, block_size, tmpbs, total, count, largest;
    block_size = tmpbs = total = count = largest = 0; 
    for (size_t i = 0; i < oshells_; i++) {
        count++;
        begin = oshell_aggs_[i];
        end = oshell_aggs_[i + 1] - 1; 
        tmpbs += end - begin + 1; 

// the name "symm" might be confused with "lr_symmetric." they are not the same
// the symm refers to whether or not the symmetry in (Q| m n) is being used
// lr_symmetric decides whether C_left_ == C_right_
        if (symm) {
            std::stringstream error;
            error << "Don't use Symm code for two-basis AOs.";
            throw PSIEXCEPTION(error.str().c_str()); 
            // in-core symmetric
            // get current cost of this block of AOs and add it to the total
            // the second buffer is accounted for with full AO_core
/*            current = symm_big_skips_[end + 1] - symm_big_skips_[begin];
            if (do_wK_) {
                current *= 3;
            }        
            total += current;*/
        } else { 
            // on-disk
            // get current cost of this block of AOs and add it to the total
            // count current twice, for both pre and post contracted buffers
            current = ob_ob_big_skips_[end + 1] - ob_ob_big_skips_[begin];
            if (do_wK_) {
                current *= 3;
            }    
            total += 2 * current;
        }            

        demand = total;
        demand += full_3index; /* bug */ /* no. not a bug. */
        demand += (hold_met_ ? naux_ * naux_ : total);
        if (demand > mem || i == oshells_ - 1) {
            if (count == 1 && i != oshells_ - 1) {
                std::stringstream error;
                error << "DFHelper: not enough memory for (o shell) AO blocking!"
                      << " required memory: " << demand * 8 / (1024 * 1024 * 1024.0) << " [GiB].";   
                throw PSIEXCEPTION(error.str().c_str());
            }
            if (demand > mem) {
                total -= current;
                tmpbs -= end - begin + 1;
                b.push_back(std::make_pair(i - count + 1, i - 1));
                i--;
            } else if (i == oshells_ - 1)
                b.push_back(std::make_pair(i - count + 1, i));

            if (largest < total) {
                largest = total;
                block_size = tmpbs;
            }
            count = 0;
            total = 0;
            tmpbs = 0;
        }
    }
    // returns pair(largest buffer size, largest block size)
    return std::make_pair(largest, block_size);
} // DFHelper_2B::fshell_blocks_for_ob_ob_AO_build

/* on variable names: b for blocks, */
/* this function serves a similar prupose to the base class function, but
   the base class function only has to work with the orbital basis

   Fundamentally, all I want to do is figure out how big I can make my temp
   arrays. What's tragic is that I have to have multiple temp arrays present
   in my program at the same time, so I have to rewrite the function to have 
   sums of all the temp sizes*/
std::tuple<size_t, size_t> DFHelper_2B::Qshell_blocks_for_JK_build(std::vector<std::pair<size_t, size_t>>& b, size_t max_nocc, bool lr_symmetric) {

    size_t T1s =  nob_ * max_nocc + ntb_ * max_nocc;
    size_t T2s = (lr_symmetric ? nob_ * nob_ + ntb_ * ntb_ + ntb_ * nob_ : nob_ * max_nocc + ntb_ * max_nocc);

    size_t T3s = std::max(nthreads_ * nob_ * nob_ , nthreads_ * nbf_ * max_nocc);
    T3s = std::max( T3s , nthreads_ * ntb_ * ntb_ );

    // 
    size_t total_AO_buffer = (AO_core_ ? big_skips_[ntb_] + ob_ao_big_skips_[nob_] + ob_ob_big_skips_[nob_] : 0 );

    size_t block_size = 0, largest = 0;
    for (size_t i = 0, tmpbs = 0, count = 1; i < Qshells_; i++, count++) {
        size_t begin  = Qshell_aggs_[i];
        size_t end = Qshell_aggs_[i + 1] - 1;

        // only really useful for the out of core method
        size_t current = (end - begin + 1) * (small_skips_[ntb_] + ob_ao_small_skips_[nob_] + ob_ob_small_skips_[nob_] );
        total_AO_buffer += (AO_core_ ? 0 : current);
        tmpbs += end - begin + 1;

        // compute total memory used by aggregate block
        size_t constraint = total_AO_buffer + T1s *tmpbs + T3s;
        constraint += (lr_symmetric ? T2s : T2s * tmpbs);
        if ( constraint > memory_ || i == Qshells_ - 1 ) {
            if (count == 1 && i == Qshells_ - 1) {
                std::stringstream error;
                error << "DFHelper: not enought memory for JK blocking";
                throw PSIEXCEPTION(error.str().c_str());
            }
            if (constraint > memory_) {
                total_AO_buffer -= (AO_core_ ? 0 : current);
                tmpbs -= end - begin + 1;
                b.push_back(std::make_pair(i - count + 1, i));
                i--;
            } else if (i == Qshells_ - 1) {
                b.push_back(std::make_pair(i - count + 1, i));
            }
            if (block_size < tmpbs) {
                largest = total_AO_buffer;
                block_size = tmpbs;
            }
            count = tmpbs = 0;
            total_AO_buffer = (AO_core_ ? big_skips_[ntb_] + ob_ao_big_skips_[nob_] + ob_ob_big_skips_[nob_] : 0 );
        }
    }
    return std::make_tuple(largest, block_size);
} // DFHelper_2B::Qshell_blocks_for_JK_build

void DFHelper_2B::prepare_AO_core() {
    timer_on("DFH_2B AO_prep");
    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    auto ao_ao_rifactory = std::make_shared<IntegralFactory>(aux_, zero, primary_, primary_);
    auto ob_ao_rifactory = std::make_shared<IntegralFactory>(aux_, zero, one_basis_, primary_);
    auto ob_ob_rifactory = std::make_shared<IntegralFactory>(aux_, zero, one_basis_, one_basis_);

    std::vector<std::shared_ptr<TwoBodyAOInt>> ao_ao_eri(nthreads_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> ob_ao_eri(nthreads_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> ob_ob_eri(nthreads_);

    std::vector<std::pair<size_t, size_t>> ao_ao_psteps;
    std::vector<std::pair<size_t, size_t>> ob_ao_fsteps;
    std::vector<std::pair<size_t, size_t>> ob_ob_fsteps;

    std::pair<size_t, size_t> ao_ao_plargest = pshell_blocks_for_AO_build(memory_, 0, ao_ao_psteps);
    std::pair<size_t, size_t> ob_ao_flargest = fshell_blocks_for_ob_ao_AO_build(memory_, 0, ob_ao_fsteps);
    std::pair<size_t, size_t> ob_ob_flargest = fshell_blocks_for_ob_ob_AO_build(memory_, 0, ob_ob_fsteps);

#pragma omp parallel num_threads(nthreads_)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        ao_ao_eri[rank] = std::shared_ptr<TwoBodyAOInt>(ao_ao_rifactory->eri());
        ob_ao_eri[rank] = std::shared_ptr<TwoBodyAOInt>(ob_ao_rifactory->eri());
        ob_ob_eri[rank] = std::shared_ptr<TwoBodyAOInt>(ob_ob_rifactory->eri());
    }

    Ppq_ = std::unique_ptr<double[]>(new double[big_skips_[nbf_]]);

//Ppq_reg_ = std::unique_ptr<double[]>(new double[big_skips_[nbf_]]);

    Pfq_ = std::unique_ptr<double[]>(new double[ob_ao_big_skips_[nob_]]);
    Pfg_ = std::unique_ptr<double[]>(new double[ob_ob_big_skips_[nob_]]);

    double* ppq = Ppq_.get();
    double* pfq = Pfq_.get();
    double* pfg = Pfg_.get();
//double* ppqr = Ppq_reg_.get();

//    std::unique_ptr<double[]> metric = std::unique_ptr<double[]>(new double[naux_ * naux_]);
    double* metp;
    mpower_ = -0.5;
//    metp = metric.get();
    metp = metric_prep_core(mpower_);

    std::unique_ptr<double[]> Qpq(new double[std::get<0>(ao_ao_plargest)]);
    double* taop = Qpq.get();
    if (do_JK_tt_ && !do_JK_oo_) {
        for (size_t i = 0; i < ao_ao_psteps.size(); i++) { 
            size_t start = std::get<0>(ao_ao_psteps[i]); 
            size_t stop = std::get<1>(ao_ao_psteps[i]);
            size_t start_func = primary_->shell(start).function_index();
            size_t stop_func = primary_->shell(stop).function_index() + primary_->shell(stop).nfunction() - 1;
            size_t begin = pshell_aggs_[start];
            size_t end = pshell_aggs_[stop + 1] - 1;
            compute_sparse_pQq_blocking_p(start, stop, taop, ao_ao_eri);
            // (P|ls) <- temp
           // C_DCOPY(big_skips_[stop_func + 1] - big_skips_[start_func], taop, 1, &ppqr[big_skips_[start_func]], 1);
            // [J^{1/2}]_{QP}(P|ls)
            contract_metric_pPq_core(start_func, stop_func, taop, metp);
        }
    }

    Qpq.reset();
    std::unique_ptr<double[]> Qfq(new double[std::get<0>(ob_ao_flargest)]);
    taop = Qfq.get();
    if (do_JK_ot_ && !do_JK_oo_) {
        for (size_t i = 0; i < ob_ao_fsteps.size(); i++) { 
            size_t start = std::get<0>(ob_ao_fsteps[i]); 
            size_t stop = std::get<1>(ob_ao_fsteps[i]);
            size_t start_func = one_basis_->shell(start).function_index();
            size_t stop_func = one_basis_->shell(stop).function_index() + one_basis_->shell(stop).nfunction() - 1;
            size_t begin = oshell_aggs_[start];
            size_t end = oshell_aggs_[stop + 1] - 1;
            compute_sparse_fPq_blocking_p(start, stop, taop, ob_ao_eri);
            // [J^{1/2}]_{QP}(P|rs)
    		contract_metric_fPq_core(start_func, stop_func, taop, metp);
        }
    }

    Qfq.reset();
    std::unique_ptr<double[]> Qfg(new double[std::get<0>(ob_ob_flargest)]);
    taop = Qfg.get();

    for (size_t i = 0; i < ob_ob_fsteps.size(); i++) { 
        size_t start = std::get<0>(ob_ob_fsteps[i]); 
        size_t stop = std::get<1>(ob_ob_fsteps[i]);
        size_t start_func = one_basis_->shell(start).function_index();
        size_t stop_func = one_basis_->shell(stop).function_index() + one_basis_->shell(stop).nfunction() - 1;
        size_t begin = oshell_aggs_[start];
        size_t end = oshell_aggs_[stop + 1] - 1;
        timer_on("DFH_2B: ERI Calculation");
        compute_sparse_fPg_blocking_p_symm(start, stop, taop, ob_ob_eri);//taop, ob_ob_eri);
        timer_off("DFH_2B: ERI Calculation");
        //contract_metric_fPq_core(start_func, stop_func, taop, metp);
        timer_on("DFH_2B: Metric Contraction");
        contract_metric_f_pPg_q_core_parsim(start_func, stop_func, taop, metp);
        timer_off("DFH_2B: Metric Contraction");
        //C_DCOPY(ob_ob_big_skips_[stop_func + 1] - ob_ob_big_skips_[start_func], taop, 1, &pfg[ob_ob_big_skips_[start_func]], 1);
        // [J^{1/2}]_{QP}(P|rq)
        //contract_metric_fPg_core( begin, end, taop, metp);
    } 
    //contract_metric_f_pPg_q_core_parsim(start_func, stop_func, taop, metp);
    Qfg.reset();
    timer_off("DFH_2B AO_prep");
} // DFHelper_2B::prepare_AO_core()

/*void DFHelper_2B::compute_dense_Pfq_blocking_p( size_t start, size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ) {

} // compute_dense_Pfq_blocking_p

void DFHelper_2B::compute_dense_Pfg_blocking_p( size_t start, size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ) {
} // compute_dense_Pfg_blocking_p */

void DFHelper_2B::compute_sparse_fPq_blocking_p( size_t start, size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ) {
    size_t begin = oshell_aggs_[start];
    size_t end = oshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    size_t start_ind = ob_ao_big_skips_[begin];

    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) {nthread = eri.size();}

    std::vector<const double*> buffer(nthread);

#pragma omp parallel num_threads(nthread)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

#pragma omp parallel for schedule(guided) num_threads(nthread)
    for (size_t PHI = start; PHI <= stop; PHI++) {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        size_t numphi = one_basis_->shell(PHI).nfunction();
        size_t indphi = one_basis_->shell(PHI).function_index();
        size_t phiind;
        for (size_t NU = 0; NU < pshells_; NU++) {
            size_t numnu = primary_->shell(NU).nfunction();
            size_t indnu = primary_->shell(NU).function_index();
            size_t nuind;
            if (!ob_ao_schwarz_shell_mask_[PHI*pshells_ + NU]) {
                continue;
            }
            for (size_t Pshell = 0; Pshell < Qshells_; Pshell++) {
                size_t Pind = aux_->shell(Pshell).function_index();
                size_t Pfuncs = aux_->shell(Pshell).nfunction();
// (AUX 0 | one primary)
                eri[rank]->compute_shell(Pshell, 0, PHI, NU);
                for ( size_t phi = 0; phi < numphi; phi++ ) {
                    phiind = indphi + phi;
                    for ( size_t nu = 0;  nu < numnu; nu++ ) {
                        nuind = indnu + nu;
                        if (ob_ao_schwarz_fun_mask_[phiind * nbf_ + nuind]) {
                            for (size_t pfunc = 0; pfunc < Pfuncs; pfunc++) {
Mp[ ob_ao_big_skips_[phiind] - start_ind + (Pind + pfunc)*ob_ao_small_skips_[phiind] +  ob_ao_schwarz_fun_mask_[phiind * nbf_ + nuind] - 1 ]
=
buffer[rank][ pfunc*numphi*numnu + phi*numnu + nu ];
                            }
                        }
                    }
                }
            }
        }// loop: size_t NU = 0; NU <= pshells_; NU++
    } // loop: size_t PHI = start; PHI <= stop; PHI++
} // compute_sparse_Pfq_blocking_p

void DFHelper_2B::compute_sparse_fPg_blocking_p( size_t start, size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ) { 
    size_t begin = oshell_aggs_[start];
    size_t end = oshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    size_t startind = ob_ob_big_skips_[begin];

    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) {nthread = eri.size();}

    std::vector<const double*> buffer(nthread);

#pragma omp parallel num_threads(nthread)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

#pragma omp parallel for num_threads(nthread) schedule(guided)
    for (size_t PHI = start; PHI <= stop; PHI++) {
        int rank = omp_get_thread_num();
        size_t indphi = one_basis_->shell(PHI).function_index(); 
        size_t numphi = one_basis_->shell(PHI).nfunction(); 
        size_t phiind;
        for (size_t GAMMA = 0; GAMMA < oshells_; GAMMA++) {
            if (!ob_ob_schwarz_shell_mask_[PHI*oshells_ + GAMMA]) {continue;}
            size_t indgamma = one_basis_->shell(GAMMA).function_index(); 
            size_t numgamma = one_basis_->shell(GAMMA).nfunction(); 
            size_t gammaind;
            for (size_t Pshell = 0; Pshell < Qshells_; Pshell++) {
                size_t indP = aux_->shell(Pshell).function_index();
                size_t numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell( Pshell, 0, PHI, GAMMA);
                for (size_t phi = 0; phi < numphi; phi++) {
                    phiind = indphi + phi;
                    for (size_t gamma = 0; gamma < numgamma; gamma++) {
                        gammaind = indgamma + gamma;
                        for (size_t pfunc = 0; pfunc < numP; pfunc++) {
                            if (ob_ob_schwarz_fun_mask_[nob_ * phiind + gammaind]) {
Mp[ob_ob_big_skips_[phiind]-startind + (indP+pfunc ) * ob_ob_small_skips_[phiind] + ob_ob_schwarz_fun_mask_[ phiind * nob_ + gammaind] - 1 ]
=
buffer[rank][ pfunc * numphi * numgamma + phi * numgamma + gamma];
                            }
                        }// size_t pfunc = 0; pfunc < numP; pfunc++
                    }// size_t gamma = 0; gamma < numgamma; gamma++
                }// size_t phi = 0; phi < numphi; phi++
            }// size_t Pshell = 0; Pshell < naux_; Pshell++
        }// size_t GAMMA = 0; GAMMA < oshells_; GAMMA++
    }// size_t PHI = start; PHI <= stop; PHI++
} // compute_sparse_Pfg_blocking_p

void DFHelper_2B::compute_sparse_fPg_blocking_p_symm( size_t start, size_t stop, double* Mp, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ) {
    size_t begin = oshell_aggs_[start];
    size_t end = oshell_aggs_[stop + 1] - 1;
    size_t block_size = end - begin + 1;
    size_t startind = ob_ao_big_skips_[one_basis_->shell(start).function_index()];

    size_t nthread = nthreads_;
    if (eri.size() != nthreads_) {nthread = eri.size();}

    std::vector<const double*> buffer(nthread);

    double* pfg = Pfg_.get();

#pragma omp parallel num_threads(nthread)
    {
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
        buffer[rank] = eri[rank]->buffer();
    }

#pragma omp parallel for num_threads(nthread) schedule(guided)
    for (size_t PHI = start; PHI <= stop; PHI++) {
        int rank = omp_get_thread_num();
        size_t indphi = one_basis_->shell(PHI).function_index(); 
        size_t numphi = one_basis_->shell(PHI).nfunction(); 
        size_t phiind;
        for (size_t GAMMA = PHI; GAMMA < oshells_; GAMMA++) {
            if (!ob_ob_schwarz_shell_mask_[PHI*oshells_ + GAMMA]) {continue;}
            size_t indgamma = one_basis_->shell(GAMMA).function_index(); 
            size_t numgamma = one_basis_->shell(GAMMA).nfunction(); 
            size_t gammaind;
            for (size_t Pshell = 0; Pshell < Qshells_; Pshell++) {
                size_t indP = aux_->shell(Pshell).function_index();
                size_t numP = aux_->shell(Pshell).nfunction();
                eri[rank]->compute_shell( Pshell, 0, PHI, GAMMA);
                for (size_t phi = 0; phi < numphi; phi++) {
                    phiind = indphi + phi;
                    for (size_t gamma = 0; gamma < numgamma; gamma++) {
                        gammaind = indgamma + gamma;
                        for (size_t pfunc = 0; pfunc < numP; pfunc++) {
                            if (ob_ob_schwarz_fun_mask_[nob_ * phiind + gammaind]) {
                                if ( ob_to_prim_[gammaind] != -1) {
Mp[

ob_ao_big_skips_[phiind] - startind + 

(indP + pfunc) * ob_ao_small_skips_[phiind] + 

ob_ao_schwarz_fun_mask_[ phiind * nbf_ + ob_to_prim_[gammaind] ] - 1  

]

=
pfg[ob_ob_big_skips_[ phiind ] + (indP+pfunc ) * ob_ob_small_skips_[ phiind ] + ob_ob_schwarz_fun_mask_[ phiind * nob_ +gammaind] - 1 ]
=
pfg[ob_ob_big_skips_[gammaind] + (indP+pfunc ) * ob_ob_small_skips_[gammaind] + ob_ob_schwarz_fun_mask_[gammaind* nob_ + phiind ] - 1 ]
=
buffer[rank][ pfunc * numphi * numgamma + phi * numgamma + gamma];
                                } else {
pfg[ob_ob_big_skips_[ phiind ] + (indP+pfunc ) * ob_ob_small_skips_[ phiind ] + ob_ob_schwarz_fun_mask_[ phiind * nob_ +gammaind] - 1 ]
=
pfg[ob_ob_big_skips_[gammaind] + (indP+pfunc ) * ob_ob_small_skips_[gammaind] + ob_ob_schwarz_fun_mask_[gammaind* nob_ + phiind ] - 1 ]
=
buffer[rank][ pfunc * numphi * numgamma + phi * numgamma + gamma];
                                }
                            }
                        }// size_t pfunc = 0; pfunc < numP; pfunc++
                    }// size_t gamma = 0; gamma < numgamma; gamma++
                }// size_t phi = 0; phi < numphi; phi++
            }// size_t Pshell = 0; Pshell < naux_; Pshell++
        }// size_t GAMMA = 0; GAMMA < oshells_; GAMMA++
    }// size_t PHI = start; PHI <= stop; PHI++
//copy the lower part of matrix into our AO_block

//size_t saboteur;
//
//for (size_t i = 0; i < nob_; i++){
//saboteur = ob_to_prim_[i];
//}
//for (size_t i = 0; i < oshells_; i++){
//saboteur = ob_to_prim_sh_[i];
//}
//for (size_t i = 0; i < nbf_; i++){
//saboteur = prim_to_ob_[i];
//}
//for (size_t i = 0; i < pshells_; i++){
//saboteur = prim_to_ob_sh_[i];
//}

//printf("png\n");
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (size_t i = one_basis_->shell(start).function_index(); 
                i < one_basis_->shell(stop).function_index() + 
                    one_basis_->shell(stop).nfunction(); i++ ) {
        for (size_t Q = 0; Q < naux_; Q++) {
            for (size_t j = 0; j <=  lastind_[i]; j++) {
                if (ob_ao_schwarz_fun_mask_[i * nbf_ + j]) {
Mp[ ob_ao_big_skips_[i] - startind + Q * 
    ob_ao_small_skips_[i] +  ob_ao_schwarz_fun_mask_[i * nbf_ +j] - 1 ]
           =
pfg[ ob_ob_big_skips_[i] + Q * ob_ob_small_skips_[i] + 
     ob_ob_schwarz_fun_mask_[ i * nob_ + prim_to_ob_[j] ] - 1 ];
                }
            }
        }
    }
} // compute_sparse_Pfg_blocking_p_symm

void DFHelper_2B::contract_metric_pPq_core(size_t func1, size_t func2, double* Qpq, double* metp) {
	double * ppq = Ppq_.get();
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (int i = func1; i <= func2; i++) {
        C_DGEMM('N', 'N', naux_, small_skips_[i], naux_, 1.0, metp, naux_, &Qpq[big_skips_[i] - big_skips_[func1] ], small_skips_[i], 0.0, &ppq[big_skips_[i]], small_skips_[i]);
    }
} // contract_metric_Ppq_core
void DFHelper_2B::contract_metric_fPq_core(size_t func1, size_t func2, double* Qfq, double* metp) {
	double * pfq = Pfq_.get();
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (int i = func1; i <= func2; i++) {
        C_DGEMM('N', 'N', naux_, ob_ao_small_skips_[i], naux_, 1.0, metp, naux_, &Qfq[ob_ao_big_skips_[i] - ob_ao_big_skips_[func1]], ob_ao_small_skips_[i], 0.0, &pfq[ob_ao_big_skips_[i]], ob_ao_small_skips_[i]);
    }
} // contract_metric_Pfq_core
void DFHelper_2B::contract_metric_fPg_core(size_t func1, size_t func2, double* Qfg, double* metp) {
 	double * pfg = Pfg_.get();
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (int i = func1; i <= func2; i++) {
        C_DGEMM('N', 'N', naux_, ob_ob_small_skips_[i], naux_, 1.0, metp, naux_, &Qfg[ob_ob_big_skips_[i]- ob_ob_big_skips_[func1]], ob_ob_small_skips_[i], 0.0, &pfg[ob_ob_big_skips_[i]], ob_ob_small_skips_[i] );
    }
} // contract_metric_Pfg_core
void DFHelper_2B::contract_metric_f_pPg_q_core_parsim(size_t start_func, size_t stop_func, double* Qfq, double* metp) {
	double * pfq = Pfq_.get();
    double * ppq = Ppq_.get();
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (int i = start_func; i <= stop_func; i++) {
        if ( ob_to_prim_[i] == -1ul ) {
            C_DGEMM('N', 'N', naux_, ob_ao_small_skips_[i], naux_, 1.0, metp, naux_, &Qfq[ob_ao_big_skips_[i] - ob_ao_big_skips_[start_func]], ob_ao_small_skips_[i], 0.0, &pfq[ob_ao_big_skips_[i]], ob_ao_small_skips_[i]);
        } else {
            C_DGEMM('N', 'N', naux_, symm_small_skips_[ ob_to_prim_[i] ], naux_, 1.0, metp, naux_, &Qfq[ob_ao_big_skips_[i] - ob_ao_big_skips_[start_func] + symm_ignored_columns_[ ob_to_prim_[i] ] ], ob_ao_small_skips_[i], 0.0, &ppq[big_skips_[ob_to_prim_[i]] + symm_ignored_columns_[ob_to_prim_[i]] ], small_skips_[ob_to_prim_[i]]);
        }
    }
#pragma omp parallel for num_threads(nthreads_) schedule(guided)
    for (size_t p = start_func; p <= stop_func ; p++) { 
        if (ob_to_prim_[p] != -1) {
            for (size_t Q = 0; Q < naux_; Q++) { 
                for (size_t q = 0; q <= ob_to_prim_[p]; q++ ) {
                    if ( ob_ao_schwarz_fun_mask_[p * nbf_ + q ] ) {
pfq[ob_ao_big_skips_[p] + Q * ob_ao_small_skips_[p] + 
    ob_ao_schwarz_fun_mask_[p * nbf_ + q] - 1 ] =
ppq[big_skips_[ob_to_prim_[p]] + Q * small_skips_[ob_to_prim_[p]] + 
    schwarz_fun_mask_[ ob_to_prim_[p] * nbf_ + q] - 1 ];
                    }
                }
                for (size_t q = ob_to_prim_[p] + 1; q < nbf_ ; q++) {
                    if ( ob_ao_schwarz_fun_mask_[p * nbf_ + q ] ) {

pfq[ob_ao_big_skips_[p] + Q * ob_ao_small_skips_[p] + 
    ob_ao_schwarz_fun_mask_[p * nbf_ + q] - 1 ]
=
pfq[ ob_ao_big_skips_[ prim_to_ob_[q] ] + Q * ob_ao_small_skips_[ prim_to_ob_[q] ] +
     ob_ao_schwarz_fun_mask_[ prim_to_ob_[q] * nbf_ + ob_to_prim_[p] ] - 1 ]
=
// I want qQ ob_to_prim[p]
ppq[ big_skips_[q] + Q * small_skips_[q] + 
     schwarz_fun_mask_[ q * nbf_ + ob_to_prim_[p] ] - 1 ]
=
ppq[ big_skips_[ob_to_prim_[p]] + Q * small_skips_[ ob_to_prim_[p] ] + 
     schwarz_fun_mask_[ob_to_prim_[p] * nbf_ + q ] - 1 ] ;

                    }
                }
            }
        }
    }
} // contract_metric_f_pPg_q_core_parsim
/* This function asks the user to supply big skips and small skips 
   which allows them (me) to transform any basis out of the AO basis!*/ 
void DFHelper_2B::first_transform_pQq(size_t bsize, size_t bcount, size_t block_size, double* Mp, double* Tp, double* Bp, size_t nbf, std::vector<size_t> sm_skips, std::vector<size_t> bg_skips, std::vector<std::vector<double>>& C_buffers, std::vector<size_t> fun_mask) {

#pragma omp parallel for schedule(guided) num_threads(nthreads_)
    for (size_t p_iter = 0; p_iter < nbf; p_iter++ ) {
        size_t sp_size = sm_skips[p_iter];
        size_t jump = big_skips_[p_iter] + bcount * sp_size; // (AO_core ? :...not implemented )
        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

/* notice we're using nbf_ and not nbf. Neither of those choices was a bug */
/* Matt didn't explicitly save the data he needed to copy over Instead, he 
   keeps track of what data is and isn't there (schwarz masks) and then he 
   goes and sees what needs to get copied every single time. It's not too 
   expensive.*/
        for (size_t q_iter = 0, sp_count = -1; q_iter < nbf_; q_iter++) {
            if (fun_mask[ p_iter * nbf + q_iter]) {
                sp_count++;
/* args: length, origin_start, origin_jump, dest_start, dest_jump*/
                C_DCOPY(bsize, Bp, 1, &C_buffers[rank][sp_count * bsize], 1);
            }
        } // size_t q_iter = 0; q_iter < nbf_; q_iter++
        //C_DGEMM( 'N', 'N', block_size, );
    }// size_t p_iter = 0; p_iter < nbf; p_iter++
}

/* Function for building all the assymmetric JK objects for now, it'll
   compute all of them by default, ...but how do we get HF to work in that
   case? I don't want to lose the sanity check */

/* another note: this is dfhelper not df_jk object. we don't actually get to 
   have J and K ourselves :( we just get to put stuff in them. (J&K aren't 
   this object's member data)
   For the user:
   This function expects, in terms of this object's member data:
   orbital coefficients Cleft & (if (!lr_symmetric) )  Cright in vectors of 
   shared matrices of dimension nbf_ by max_nocc 
   SCF density matrix D of dimension nbf_ by nbf_
   AO        by AO        coulomb matrix J_pq of dimension nbf_ by nbf_ 
   one basis by AO        coulomb matrix J_fq of dimension nob_ by nbf_ 
   one basis by one basis coulomb matrix J_fg of dimension nob_ by nob_ 
   AO        by AO        exchange matrix K_pq of dimension nbf_ by nbf_ 
   one basis by AO        exchange matrix K_fq of dimension nob_ by nbf_ 
   one basis by one basis exchange matrix K_fg of dimension nob_ by nob_  */
void DFHelper_2B::build_2B_JK_naive(std::vector<SharedMatrix> Cleft,
                              std::vector<SharedMatrix> Cright,
                              std::vector<SharedMatrix> D,
                              std::vector<SharedMatrix> J_pq,
                              std::vector<SharedMatrix> J_fq,
                              std::vector<SharedMatrix> J_fg,
                              std::vector<SharedMatrix> K_pq,
                              std::vector<SharedMatrix> K_fq,
                              std::vector<SharedMatrix> K_fg,
                              size_t max_nocc, bool do_J, bool do_K, 
                              bool lr_symmetric) {

    /* dummy vector to use matt's code*/ 
    /* std::vector<SharedMatrix> wk; */

    /* wk.push_back(Matrix()); */

    /* we need the AO integrals (P|pq)*/
    double* ppq = Ppq_.get();
    /* we need the AO integrals (P|fq)*/
    double* pfq = Pfq_.get();
    /* we need the AO integrals (P|fg)*/
    double* pfg = Pfg_.get();


//double* ppqr = Ppq_reg_.get();

//double* met_inv = metric_inverse_prep_core();

double* metp = metric_prep_core(-0.5);
double* met_lu_p = Metric_LU_.get();
int* pert = MET_LU_PERT_.get();


    std::unique_ptr<double[]> Q_vec(new double[naux_]); 
    double* qvec = Q_vec.get();

std::unique_ptr<double[]> Q_vec2(new double[naux_]); 
double* qvec2 = Q_vec2.get();

    std::unique_ptr<double[]> T1_pao(new double[max_nocc*naux_*nbf_]);
    std::unique_ptr<double[]> T1_fao(new double[max_nocc*naux_*nob_]);

    std::unique_ptr<double[]> T2_pao;
    std::unique_ptr<double[]> T2_fao;

    double* t1_pp = T1_pao.get();
    double* t1_fp = T1_fao.get();
    double* t2_pp;
    double* t2_fp;
    if (!lr_symmetric) {
        /* max_nocc is correct... enough, not the right amount*/
        T2_pao = std::make_unique<double[]>(max_nocc*naux_*nbf_);
        T2_fao = std::make_unique<double[]>(max_nocc*naux_*nob_);

        t2_pp = T2_pao.get();
        t2_fp = T2_fao.get();
    } else {
        t2_pp = T1_pao.get();
        t2_fp = T1_fao.get();
    }

    size_t l_nocc;
    size_t r_nocc;

    double* dp;//  = D[0]->pointer();
    double* clp;// = Cleft
    double* crp;
    double* j_pq;
    double* j_fq;
    double* j_fg;
    double* k_pq;
    double* k_fq;
    double* k_fg;

    /* why instantiate a bunch of zeros? we're going to deep copy out of 
     * this vector. A smart programmer will try to copy a whole cache-line, so
     * I'm including a little extra space in case they got a little sloppy copying
     * the cache-line */
    double dub_zero[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };


    for ( int N = 0; N <D.size(); N++ ) {
        dp = D[N]->pointer()[0];
        clp = Cleft[N]->pointer()[0];
        crp = Cright[N]->pointer()[0];

        l_nocc = Cleft[N]->coldim(0);
        r_nocc = Cright[N]->coldim(0);

        j_pq = J_pq[N]->pointer()[0];
        j_fq = J_fq[N]->pointer()[0];
        j_fg = J_fg[N]->pointer()[0];
        k_pq = K_pq[N]->pointer()[0];
        k_fq = K_fq[N]->pointer()[0];
        k_fg = K_fg[N]->pointer()[0];

        /* There's a bit of an issue here... we need the contracted pfq 
           the same time as the contracted ppq. Just a bit of an issue, we can
           move across memory (kind of ) half as quickly */
    
       /* we're going to build both the J's in this block. Don't get too upset
          this is going to be replaced by some cleaner code soon enough. 
    
          I'm going to write my own function for the J construction. We get to
          reuse the precursor BD to form J I think it's lovely.*/
    
       /* from the expression BBD, the first conctraction v <- BD*/
       /* There're a million things wrong with this line. Among them, that in
          an actual sparse case, we'd have to loop over all basis functions
          and prune down the density matrix... but it's worth it for the 
          savings in exchange*/
       /* Save a tree... Use Density Fitting! we got to reuse qvec
          even at this juncture, though, we still have to make a lot of
          seperat gemv's ... each in its own loop. drats!*/


       /* number of one basis funcs. nob_ number two basis "" ntb_ */

       /*  j_pq  over basis function iter*/
       /* fill with 0's :) */
        C_DCOPY(naux_, dub_zero, 0, qvec, 1);
        /* prepare the contraction V_Q <- (Q|mn)D_mn */
//        for (size_t nbf_iter = 0; nbf_iter < nbf_; nbf_iter++ ) {
//            C_DGEMV( 'N', naux_, nbf_, 1.0, ppq + nbf_iter * nbf_ * naux_, nbf_, dp + nbf_iter * nbf_, 1, 1.0, qvec, 1);
//        }


        // calculates V <- (P|ls)D_{ls}
        for (size_t nbf_iter = 0; nbf_iter < nbf_; nbf_iter++ ) {
            C_DGEMV( 'N', naux_, nbf_, 1.0, ppq + nbf_iter * nbf_ * naux_, nbf_, dp + nbf_iter * nbf_, 1, 1.0, qvec, 1);
        }

//ppqr

        // Metric Contraction for the Coulomb Matrix
        //C_DGETRS('N', naux_, 1, met_lu_p, naux_, pert, qvec, naux_);
        // somthing unholy...
        C_DGEMV('N', naux_, naux_, 1.0, metp, naux_, qvec, 1, 0.0, qvec2, 1);
//        C_DCOPY(naux_, qvec2, 1, qvec, 1);

/*        for (int bf_iter = 0; bf_iter < nbf_; bf_iter++ ) {
            C_DGEMV('T', naux_, nbf_, 1.0, ppqr + naux_ * nbf_ * bf_iter, nbf_, qvec2, 1, 0.0, j_pq + nbf_ * bf_iter , 1);
        } */


        for (int of_iter = 0; of_iter < nob_; of_iter++ ) {
            C_DGEMV('T', naux_, nbf_, 1.0, pfq + naux_ * nbf_ * of_iter, nbf_, qvec, 1, 0.0, j_fq + nbf_ * of_iter, 1 );
        }


        for (int of_iter = 0; of_iter < nob_; of_iter++ ) {
            C_DGEMV('T', naux_, nob_, 1.0, pfg + naux_ * nob_ * of_iter, nob_, qvec2, 1, 0.0, j_fg + nob_ * of_iter, 1 );
        }

        /* These two DGEMM's are the "first contraction" i.e. (Q|ma) <- (Q|mn)C_la
           We're making two contractions... for THREE exchange matrices!!! how?   
           well... we get to use (Q|fa) for K_fg and for K_fq  
                 & we get to use (Q|pa) for K_fq and for K_pq                   */

        C_DGEMM('N', 'N', nbf_ * naux_, l_nocc, nbf_, 1.0, ppq, nbf_, clp, l_nocc, 0.0, t1_pp, l_nocc);
        C_DGEMM('N', 'N', nob_ * naux_, l_nocc, nbf_, 1.0, pfq, nbf_, clp, l_nocc, 0.0, t1_fp, l_nocc);

 
        if (!lr_symmetric) {
            C_DGEMM('N', 'N', nbf_ * naux_, r_nocc, nbf_, 1.0, ppq, nbf_, crp, r_nocc, 0.0, t2_pp, r_nocc);
            C_DGEMM('N', 'N', nob_ * naux_, r_nocc, nbf_, 1.0, pfq, nbf_, crp, r_nocc, 0.0, t2_fp, r_nocc);
        } 

        /* These three contractions form the three exchange matrices we've promised
           the user. in order K_pq, K_fq, K_fg */
        C_DGEMM('N', 'T', nbf_, nbf_, naux_ * l_nocc, 1.0, t1_pp, naux_ * l_nocc, t2_pp, naux_ * r_nocc, 0.0, k_pq, nbf_);
        C_DGEMM('N', 'T', nob_, nbf_, naux_ * l_nocc, 1.0, t1_fp, naux_ * l_nocc, t2_pp, naux_ * r_nocc, 0.0, k_fq, nbf_);
        C_DGEMM('N', 'T', nob_, nob_, naux_ * l_nocc, 1.0, t1_fp, naux_ * l_nocc, t2_fp, naux_ * r_nocc, 0.0, k_fg, nob_);
    }
    
}

void DFHelper_2B::build_2B_JK(std::vector<SharedMatrix> Cleft,
                              std::vector<SharedMatrix> Cright,
                              std::vector<SharedMatrix> D,
                              std::vector<SharedMatrix> J_pq,
                              std::vector<SharedMatrix> J_fq,
                              std::vector<SharedMatrix> J_fg,
                              std::vector<SharedMatrix> K_pq,
                              std::vector<SharedMatrix> K_fq,
                              std::vector<SharedMatrix> K_fg,
                              size_t max_nocc, bool do_J, bool do_K, 
                              bool lr_symmetric) {

    std::vector<std::pair<size_t, size_t>> Qsteps;
    std::tuple<size_t, size_t> info = Qshell_blocks_for_JK_build(Qsteps, max_nocc, lr_symmetric);
    size_t tots = std::get<0>(info);
    size_t totsb = std::get<1>(info);

    /* dummy vector to use matt's code*/ 
    /* std::vector<SharedMatrix> wk; */

    /* wk.push_back(Matrix()); */

    /* we need the AO integrals (P|pq)*/
    double* ppq = Ppq_.get();
    /* we need the AO integrals (P|fq)*/
    double* pfq = Pfq_.get();
    /* we need the AO integrals (P|fg)*/
    double* pfg = Pfg_.get();


    //double* ppqr = Ppq_reg_.get();

//double* met_inv = metric_inverse_prep_core();

    double* metp = metric_prep_core(-0.5);
    double* met_lu_p = Metric_LU_.get();
    int* pert = MET_LU_PERT_.get();

//    std::unique_ptr<double[]> C_Buffers(new double[nthreads_ * nbf_ * nbf_]);
//    double* c_buffers = C_Buffers.get();

    std::vector<std::vector<double>> c_buffers;
    c_buffers.resize(nthreads_);

#pragma omp parallel
{
    int rank = 0;
#ifdef _OPENMP
    rank = omp_get_thread_num();
#endif
    c_buffers[rank].resize(nbf_ * nbf_ );
}

    std::unique_ptr<double[]> Q_vec(new double[naux_ * nthreads_]);
    double* qvec  = Q_vec.get();
    std::unique_ptr<double[]> Q_vec2(new double[naux_ * nthreads_]);
    double* qvec2 = Q_vec2.get();

    size_t T1_mem = (!max_nocc ? totsb : totsb * max_nocc);
    T1_mem = std::max(T1_mem * nob_, naux_ * nthreads_);
    T1_mem = std::max(T1_mem, nob_ * nob_);
    T1_mem = std::max(T1_mem, ntb_ * ntb_);

    size_t T2_mem = (!max_nocc ? totsb : totsb * max_nocc);
    T2_mem = std::max(T2_mem * ntb_, ntb_ * ntb_);
    T2_mem = std::max(T2_mem, totsb);

    std::unique_ptr<double[]> T1_pao(new double[T2_mem]);
    std::unique_ptr<double[]> T1_fao(new double[T1_mem]);

    std::unique_ptr<double[]> T2_pao;
    std::unique_ptr<double[]> T2_fao;


    double* t1_pp = T1_pao.get();
    double* t1_fp = T1_fao.get();
    double* t2_pp;
    double* t2_fp;
    if (!lr_symmetric) {
        /* max_nocc is correct... enough, not the right amount*/
        T2_pao = std::make_unique<double[]>(T2_mem);
        T2_fao = std::make_unique<double[]>(T1_mem);

        t2_pp = T2_pao.get();
        t2_fp = T2_fao.get();
    } else {
        t2_pp = T1_pao.get();
        t2_fp = T1_fao.get();
    }

    size_t l_nocc;
    size_t r_nocc;

//    double* dp;//  = D[0]->pointer();
    double* clp;// = Cleft
    double* crp;
    double* k_pq;
    double* k_fq;
    double* k_fg;

    /* why instantiate a bunch of zeros? we're going to deep copy out of 
     * this vector. A smart programmer will try to copy a whole cache-line, so
     * I'm including a little extra space in case they got a little sloppy copying
     * the cache-line */
    double dub_zero[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    /* this function can't block due to the working equation we've chosen 
       for the coulomb matrix. The last two arguments are superfluous*/
    timer_on("DFH_2B J");
    if (AO_core_) {
        compute_J_2B_core(D, J_pq, J_fq, J_fg, c_buffers, ppq, pfq, pfg, t1_fp, qvec, qvec2, metp, 0, naux_);
    }
    timer_off("DFH_2B J");

    for (size_t bl_iter = 0, bcount=0; bl_iter < Qsteps.size(); bl_iter++) {
        size_t start = std::get<0>(Qsteps[bl_iter]);
        size_t stop = std::get<1>(Qsteps[bl_iter]);
        size_t begin = Qshell_aggs_[start];
        size_t end = Qshell_aggs_[stop + 1] - 1;
        size_t block_size = end - begin + 1;

        compute_K_2B(Cleft, Cright, K_pq, K_fq, K_fg, t1_pp, t2_pp, t1_fp, t2_fp, ppq, pfq, bcount, block_size, lr_symmetric, c_buffers);

        bcount += block_size;
    }




//    for ( int N = 0; N <D.size(); N++ ) {
/*        dp = D[N]->pointer()[0];
        clp = Cleft[N]->pointer()[0];
        crp = Cright[N]->pointer()[0];

        l_nocc = Cleft[N]->coldim(0);
        r_nocc = Cright[N]->coldim(0);

        j_pq = J_pq[N]->pointer()[0];
        j_fq = J_fq[N]->pointer()[0];
        j_fg = J_fg[N]->pointer()[0];
        k_pq = K_pq[N]->pointer()[0];
        k_fq = K_fq[N]->pointer()[0];
        k_fg = K_fg[N]->pointer()[0]; */

        /* There's a bit of an issue here... we need the contracted pfq 
           the same time as the contracted ppq. Just a bit of an issue, we can
           move across memory (kind of ) half as quickly */
    
       /* we're going to build both the J's in this block. Don't get too upset
          this is going to be replaced by some cleaner code soon enough. 
    
          I'm going to write my own function for the J construction. We get to
          reuse the precursor BD to form J I think it's lovely.*/
    
       /* from the expression BBD, the first conctraction v <- BD*/
       /* There're a million things wrong with this line. Among them, that in
          an actual sparse case, we'd have to loop over all basis functions
          and prune down the density matrix... but it's worth it for the 
          savings in exchange*/
       /* Save a tree... Use Density Fitting! we got to reuse qvec
          even at this juncture, though, we still have to make a lot of
          seperat gemv's ... each in its own loop. drats!*/


       /* number of one basis funcs. nob_ number two basis "" ntb_ */

       /*  j_pq  over basis function iter*/
       /* fill with 0's :) */
//        C_DCOPY(naux_, dub_zero, 0, qvec, 1);
        /* prepare the contraction V_Q <- (Q|mn)D_mn */
//        for (size_t nbf_iter = 0; nbf_iter < nbf_; nbf_iter++ ) {
//            C_DGEMV( 'N', naux_, nbf_, 1.0, ppq + nbf_iter * nbf_ * naux_, nbf_, dp + nbf_iter * nbf_, 1, 1.0, qvec, 1);
//        }


        // calculates V <- (P|ls)D_{ls}
//        for (size_t nbf_iter = 0; nbf_iter < nbf_; nbf_iter++ ) {
//            C_DGEMV( 'N', naux_, nbf_, 1.0, ppqr + nbf_iter * nbf_ * naux_, nbf_, dp + nbf_iter * nbf_, 1, 1.0, qvec, 1);
//        }

//ppqr

        // Metric Contraction for the Coulomb Matrix
//        C_DGETRS('N', naux_, 1, met_lu_p, naux_, pert, qvec, naux_);
        // somthing unholy...
//        C_DGEMV('N', naux_, naux_, 1.0, metp, naux_, qvec, 1, 0.0, qvec2, 1);
//        C_DCOPY(naux_, qvec2, 1, qvec, 1);

//        for (int bf_iter = 0; bf_iter < nbf_; bf_iter++ ) {
//            C_DGEMV('T', naux_, nbf_, 1.0, ppqr + naux_ * nbf_ * bf_iter, nbf_, qvec, 1, 0.0, j_pq + nbf_ * bf_iter , 1);
//        }
//
//
//        for (int of_iter = 0; of_iter < nob_; of_iter++ ) {
//            C_DGEMV('T', naux_, nbf_, 1.0, pfq + naux_ * nbf_ * of_iter, nbf_, qvec, 1, 0.0, j_fq + nbf_ * of_iter, 1 );
//        }
//
//
//        for (int of_iter = 0; of_iter < nob_; of_iter++ ) {
//            C_DGEMV('T', naux_, nob_, 1.0, pfg + naux_ * nob_ * of_iter, nob_, qvec, 1, 0.0, j_fg + nob_ * of_iter, 1 );
//        }

        /* These two DGEMM's are the "first contraction" i.e. (Q|ma) <- (Q|mn)C_la
           We're making two contractions... for THREE exchange matrices!!! how?   
           well... we get to use (Q|fa) for K_fg and for K_fq  
                 & we get to use (Q|pa) for K_fq and for K_pq                   */

/*        C_DGEMM('N', 'N', nbf_ * naux_, l_nocc, nbf_, 1.0, ppq, nbf_, clp, l_nocc, 0.0, t1_pp, l_nocc);
        C_DGEMM('N', 'N', nob_ * naux_, l_nocc, nbf_, 1.0, pfq, nbf_, clp, l_nocc, 0.0, t1_fp, l_nocc);

 
        if (!lr_symmetric) {
            C_DGEMM('N', 'N', nbf_ * naux_, r_nocc, nbf_, 1.0, ppq, nbf_, crp, r_nocc, 0.0, t2_pp, r_nocc);
            C_DGEMM('N', 'N', nob_ * naux_, r_nocc, nbf_, 1.0, pfq, nbf_, crp, r_nocc, 0.0, t2_fp, r_nocc);
        }   */

        /* These three contractions form the three exchange matrices we've promised
           the user. in order K_pq, K_fq, K_fg */
        /* C_DGEMM('N', 'T', nbf_, nbf_, naux_ * l_nocc, 1.0, t1_pp, naux_ * l_nocc, t2_pp, naux_ * r_nocc, 0.0, k_pq, nbf_);
        C_DGEMM('N', 'T', nob_, nbf_, naux_ * l_nocc, 1.0, t1_fp, naux_ * l_nocc, t2_pp, naux_ * r_nocc, 0.0, k_fq, nbf_);
        C_DGEMM('N', 'T', nob_, nob_, naux_ * l_nocc, 1.0, t1_fp, naux_ * l_nocc, t2_fp, naux_ * r_nocc, 0.0, k_fg, nob_); */
//    }
    
}
/* A note to the user about the many arguments to this function:
   This function is meant to be used in perturbation theories, so we don't
   deal with a convergent method. Instead, we get one set of orbitals out of
   which to form Coulomb and Exchange operators. 
   D is the converved SCF density matrix
   J_pq is the coulomb matrix in the two basis.
   J_fq is the coulomb matrix one_basis X two_basis.
   J_fg is the coulomb matrix in the one basis
   */
void DFHelper_2B::compute_J_2B_core(std::vector<SharedMatrix> D, 
                                    std::vector<SharedMatrix> J_pq,
                                    std::vector<SharedMatrix> J_fq,
                                    std::vector<SharedMatrix> J_fg,
                                    std::vector<std::vector<double>>& dbuffers,
                                    double* Qpq,
                                    double* Qfq,
                                    double* Pfg,
                                    double* T1p,
                                    double* qvec,
                                    double* qvec2,
                                    double* metp,
                                    size_t bcount,
                                    size_t block_size){
//JSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOBJSOB
int n = 0;
    double ZERO[16] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double* zerop = &ZERO[0];
// does is make sense for all the variables compute_J_2B to be computed
// inside the loop? yes. They should have their own AO spaces
    for ( size_t N = 0; N < D.size(); N++ ){
n = 1;
        double* Dp = D[N]->pointer()[0];
        double* j_pq = J_pq[N]->pointer()[0];
        double* j_fq = J_fq[N]->pointer()[0];
        double* j_fg = J_fg[N]->pointer()[0];
        C_DCOPY(block_size, zerop, 0, qvec, 1);
        C_DCOPY(block_size * nthreads_, zerop, 0, T1p, 1);
// V = (Q|mn) D[N]_{mn}
// V = B_m^Q_n D[N]_{m n}
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for ( size_t bf_iter = 0; bf_iter < nbf_; bf_iter++ ){
            size_t pq_size = small_skips_[bf_iter];
            size_t jump = (AO_core_ ? big_skips_[bf_iter] + bcount * pq_size : (big_skips_[bf_iter] * block_size) / naux_ );
            int rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif
            for (size_t bf2 = 0, count = -1; bf2 < nbf_; bf2++) {
                if (schwarz_fun_mask_[ bf_iter * nbf_ + bf2 ]) {
                    count++;
                    dbuffers[rank][count] = Dp[nbf_ * bf_iter + bf2];
                }
            }

            C_DGEMV('N', block_size, pq_size, 1.0, &Qpq[jump], pq_size, &dbuffers[rank][0], 1, 1.0, &T1p[rank * block_size], 1);
        }
        for (size_t i = 0; i < nthreads_; i++) {
            for (size_t bf_iter = 0; bf_iter < block_size; bf_iter++) {
                qvec[bf_iter] += T1p[ i * block_size + bf_iter];
            }
        }

// V2 = (P|Q) V_Q
        C_DGEMV('N', block_size, naux_, 1.0, metp + naux_ * bcount, naux_, qvec, 1, 0.0, qvec2, 1); 


// Jpq = B_pQq V_Q or Jfq = V_Q or J_fg = A_fPg V2_P
        C_DCOPY(nbf_ * nbf_, zerop, 0, T1p, 1);
        C_DCOPY(nbf_ * nbf_, zerop, 0, j_pq,1);

#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t bf_iter = 0; bf_iter < nbf_; bf_iter++) {
            int rank = 0;
#ifdef _OPENMP
            rank = omp_get_thread_num();
#endif
            size_t sp_size = small_skips_[bf_iter];
            size_t jump = (AO_core_ ? big_skips_[bf_iter] + bcount * sp_size : (big_skips_[bf_iter] * block_size)/naux_ );
            C_DGEMV('T', block_size, sp_size, 1.0, Qpq + jump, sp_size, qvec, 1, 0.0, &T1p[ nbf_ * bf_iter ], 1);
        }

        for (int bf_iter = 0; bf_iter < nbf_; bf_iter++) {
            for (int bf2 = 0, count = -1; bf2 < nbf_; bf2++) {
                if (schwarz_fun_mask_[ bf_iter * nbf_ + bf2]) {
                    count++;
                    j_pq[bf_iter * nbf_ + bf2] += T1p[nbf_ * bf_iter + count];
                }
            }
        }

        C_DCOPY(nob_ * nbf_, zerop, 0, T1p, 1);
        C_DCOPY(nob_ * nbf_, zerop, 0, j_fq,1);
        if (do_JK_ot_ && !do_JK_oo_) {
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
            for (size_t ob_iter = 0; ob_iter < nob_; ob_iter++) {
                int rank = 0;
#ifdef _OPENMP
                rank = omp_get_thread_num();
#endif
                size_t sp_size = ob_ao_small_skips_[ob_iter];
                size_t jump = (AO_core_ ? ob_ao_big_skips_[ob_iter] + bcount * sp_size : (ob_ao_big_skips_[ob_iter] * block_size) / naux_ );
                C_DGEMV('T', block_size, sp_size, 1.0, Qfq + jump, sp_size, qvec, 1, 0.0, &T1p[nbf_ * ob_iter], 1);
            }

            for (int ob_iter = 0; ob_iter < nob_; ob_iter++) {
                for (int bf2 = 0, count = -1; bf2 < nbf_; bf2++) {
                    if (ob_ao_schwarz_fun_mask_[ ob_iter * nbf_ + bf2]) {
                        count++;
                        j_fq[ob_iter * nbf_ + bf2] += T1p[ob_iter * nbf_ + count];
                    }
                }
            }
        }

        C_DCOPY(nob_ * nob_ , zerop, 0, T1p, 1);
        C_DCOPY(nob_ * nob_ , zerop, 0, j_fg, 1);
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
        for (size_t ob_iter = 0; ob_iter < nob_; ob_iter++) {
             int rank = 0;
#ifdef _OPENMP
             rank = omp_get_thread_num();
#endif
             size_t sp_size = ob_ob_small_skips_[ob_iter];
             size_t jump = (AO_core_ ? ob_ob_big_skips_[ob_iter] + bcount * sp_size : (ob_ob_big_skips_[ob_iter] * block_size) / naux_ );
             C_DGEMV('T', block_size, sp_size, 1.0, Pfg + jump, sp_size, qvec2, 1, 0.0, &T1p[nob_ * ob_iter], 1);
        }

        for (int ob_iter = 0; ob_iter < nob_; ob_iter++) {
            for (int bf2 = 0, count = -1; bf2 < nob_; bf2++) {
                if (ob_ob_schwarz_fun_mask_[ ob_iter * nob_ + bf2]) {
                    count++;
                    j_fg[ob_iter * nob_ + bf2] += T1p[ob_iter * nob_ + count];
                }
            }
        } 

        if (do_JK_oo_ && do_JK_ot_) {
            for (int i = 0; i < nob_; i++ ) {
                for (int j = 0; j < nbf_; j++) {
                    j_fq[i * nbf_ + j ] = j_fg[ i * nob_ + prim_to_ob_[j] ];
                }
            }
        }
    }
} // DFHelper_2B::compute_J_2B()

void DFHelper_2B::compute_K_2B(std::vector<SharedMatrix> Cleft,
                               std::vector<SharedMatrix> Cright,
                               std::vector<SharedMatrix> K_pq,
                               std::vector<SharedMatrix> K_fq,
                               std::vector<SharedMatrix> K_fg,
                               double* T1t_p,
                               double* T2t_p,
                               double* T1o_p,
                               double* T2o_p,
                               double* Mtp,
                               double* Mop,
                               size_t bcount,
                               size_t block_size,
                               bool lr_symmetric,
                               std::vector<std::vector<double>>& C_buffers) {

    for (size_t i = 0; i < Cleft.size(); i++) {
        size_t nocc = Cleft[i]->colspi()[0];
        if (!nocc) {
            continue;
        }

        double* Clp = Cleft[i]->pointer()[0];
        double* Crp = Cright[i]->pointer()[0];
        double* Kpqp = K_pq[i]->pointer()[0];
        double* Kfqp = K_fq[i]->pointer()[0];
        double* Kfgp = K_fg[i]->pointer()[0];

        timer_on("DFH_2B K1");
        if (!do_JK_oo_ && (do_JK_tt_ || do_JK_ot_) ) {    
            first_transform_pQq(nocc, bcount, block_size, Mtp, T1t_p, Clp, ntb_, nbf_, big_skips_, small_skips_, schwarz_fun_mask_, C_buffers);
        }
        if (do_JK_oo_) {
            first_transform_pQq(nocc, bcount, block_size, Mop, T1o_p, Clp, nob_, nbf_, ob_ao_big_skips_, ob_ao_small_skips_, ob_ao_schwarz_fun_mask_, C_buffers);
        }

        if (!lr_symmetric) {
            if (!do_JK_oo_ && (do_JK_tt_ || do_JK_ot_) ) {    
                first_transform_pQq(nocc, bcount, block_size, Mtp, T2t_p, Crp, ntb_, nbf_, big_skips_, small_skips_, schwarz_fun_mask_, C_buffers);
                }
            if (do_JK_oo_) {
                first_transform_pQq(nocc, bcount, block_size, Mop, T2o_p, Crp, nob_, nbf_, ob_ao_big_skips_, ob_ao_small_skips_, ob_ao_schwarz_fun_mask_, C_buffers);
            }
        }
        timer_off("DFH_2B K1");

        timer_on("DFH_2B K2");
        if (do_JK_tt_ && !do_JK_oo_) {
            C_DGEMM('N', 'T', ntb_, ntb_, block_size * nocc, 1.0, T1t_p, nocc * block_size, T2t_p, nocc * block_size, 1.0, Kpqp, ntb_);
        }

        if ( do_JK_ot_ && !do_JK_oo_ ) {
            C_DGEMM('N', 'T', nob_, ntb_, block_size * nocc, 1.0, T1o_p, nocc * block_size, T2t_p, nocc * block_size, 1.0, Kfqp, ntb_);
        }

        if (do_JK_oo_) {
            C_DGEMM('N', 'T', nob_, nob_, block_size * nocc, 1.0, T1o_p, nocc * block_size, T2o_p, nocc * block_size, 1.0, Kfgp, nob_);
        }
        timer_off("DFH_2B K2");


        if (do_JK_oo_ && do_JK_ot_) {
            for (int i = 0; i < nob_; i++ ) {
                for (int j = 0; j < nbf_; j++) {
                    Kfqp[i * nbf_ + j ] = Kfgp[ i * nob_ + prim_to_ob_[j] ];
                }
            }
        }

    }
    
} //DFHelper_2B::compute_K_2B
// look at the arguments to the dgemm if you need to know what the arguments to
//   first_transform mean
void DFHelper_2B::first_transform_pQq(size_t bsize, 
                                      size_t bcount, 
                                      size_t block_size, 
                                      double* Mp, 
                                      double* Tp, 
                                      double* Bp, 
                                      size_t rnbf, 
                                      size_t cnbf, 
                                      std::vector<size_t> big_skips, 
                                      std::vector<size_t> small_skips, 
                                      std::vector<size_t> schwarz_fun_mask, 
                                      std::vector<std::vector<double>>& C_buffers) {
#pragma omp parallel for schedule(guided) num_threads(nthreads_)
    for (size_t k = 0; k < rnbf; k++) {
        size_t sp_size = small_skips[k];
        size_t jump = (AO_core_ ? big_skips[k] + bcount * sp_size : (big_skips[k] * block_size) / naux_);

        int rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif
// we're just taking rows out of the coefficient matrix
        for (size_t m = 0, sp_count = -1; m < cnbf; m++) {
            if (schwarz_fun_mask[ k * cnbf + m ]) {
                sp_count++;
/* sp_count will go up to sp_size -1, but we don't know when it does at compile
      time, so we need the conditional increment instead */

                C_DCOPY(bsize, &Bp[m * bsize], 1, &C_buffers[rank][sp_count * bsize], 1);
            }
        }
// perform the actual transformation
        C_DGEMM('N', 'N', block_size, bsize, sp_size, 1.0, &Mp[jump], sp_size, &C_buffers[rank][0], bsize, 0.0, &Tp[k * block_size * bsize], bsize );
    }
} // DFHelper_2B::first_transform_pQq

void DFHelper_2B::set_do_JK_tt(bool do_JK_tt) {do_JK_tt_ = do_JK_tt;}
bool DFHelper_2B::get_do_JK_tt() {return do_JK_tt_;} 
void DFHelper_2B::set_do_JK_ot(bool do_JK_ot) {do_JK_ot_ = do_JK_ot;} 
bool DFHelper_2B::get_do_JK_ot() {return do_JK_ot_;} 
void DFHelper_2B::set_do_JK_oo(bool do_JK_oo) {do_JK_oo_ = do_JK_oo;}
bool DFHelper_2B::get_do_JK_oo() {return do_JK_oo_;} 

} // namespace psi

