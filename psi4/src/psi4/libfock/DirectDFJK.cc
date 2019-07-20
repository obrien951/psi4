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

#include "psi4/libfock/jk.h"

#include "psi4/lib3index/3index.h"
#include "psi4/lib3index/dftensor.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/matrix.h"
#include "psi4/pragma.h"
#include "psi4/libqt/qt.h"

#include <cstdlib>
#include <memory>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/aiohandler.h"
#include <sstream>
#include "psi4/libpsi4util/PsiOutStream.h"



namespace psi {

DirectDFJK::DirectDFJK(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary )
    : JK(primary), auxiliary_(auxiliary)  {
// Not much happens in common init, just that the number of elements in the
//  basis sets gets saved as private member data in the DirectDFJK object.
//  However, the JK() constructor gets called as well. Not to fear!
//  In our case, all that does is set a couple booleans. 
    common_init();
}

DirectDFJK::~DirectDFJK() {}

//It's hard to say what does and doesn't belong in this function. I am
// open to suggestions as to what should change. JSOB
void DirectDFJK::common_init(){
    nbf_ = primary_->nbf();
    naux_ = auxiliary_->nbf();

	p_shells_ = primary_->nshell();
	Q_shells_ = auxiliary_->nshell();
}

//Used as testing machinery. It's nice to have a public function
//to be able to look at private member data from pyside.
//I plan on nuking this as soon as I get this project done,
//don't you worry. JSOB

void DirectDFJK::pytemplate(){
	int* five = (int*) malloc(sizeof(int));
	int* twenty_four = (int*) malloc(sizeof(int));
                                          
	five[0] = 5;
	twenty_four[0] = 24;

//	double* jjp;
//	double* jkp;
//	double* mjp;
//	double* mkp;
//	double* jdp;
//	double* kdp;

//	SharedMatrix Joe__J = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix Joe__K = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix Matt_J = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix Matt_K = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix jdiff = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix kdiff = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);

//	jjp = Joe__J->pointer()[0];
//	jkp = Joe__K->pointer()[0];
//	mjp = Matt_J->pointer()[0];
//	mkp = Matt_K->pointer()[0];
//	jdp = jdiff->pointer()[0];
//	kdp = kdiff->pointer()[0];

//	Matt_J->load("/theoryfs2/ds/obrien/Debug/Psi4/memdfjk_j.txt");
//	Matt_K->load("/theoryfs2/ds/obrien/Debug/Psi4/memdfjk_k.txt");

	SharedMatrix tmp_c = std::make_shared<Matrix>(1, twenty_four, five, 0);
	
	tmp_c->load("/theoryfs2/ds/obrien/Debug/Psi4/C_memdfjk.txt");

	C_left_ao_.push_back(tmp_c);
	C_right_ao_.push_back(tmp_c);
	compute_D_ao(0);

	std::stringstream J_name;
	std::stringstream K_name;
   	J_name << "J " << 0 << " (AO)";
   	K_name << "K " << 0 << " (AO)";
	J_ao_.push_back(std::make_shared<Matrix>(J_name.str(), (int) nbf_, (int) nbf_));
	K_ao_.push_back(std::make_shared<Matrix>(K_name.str(), (int) nbf_, (int) nbf_));
	
	num_blocks_ = 2;

	Block_funcs_.resize(0);

	Block_funcs_.push_back(0);
	Block_funcs_.push_back(0);

	Shell_starts_.push_back(0);
	Shell_starts_.push_back((Q_shells_/2) + 1);

	Shell_stops_.push_back(Q_shells_/2);
	Shell_stops_.push_back(Q_shells_ - 1);

	for (size_t i = 0; i <= Q_shells_/2; i++) {
		Block_funcs_[0] += auxiliary_->shell(i).nfunction();
	}

	for (size_t i = (Q_shells_/2) + 1; i < Q_shells_; i++) {
		Block_funcs_[1] += auxiliary_->shell(i).nfunction();
	}

	biggest_block_ = std::max( Block_funcs_[1], Block_funcs_[0]);

	biggest_block_ = biggest_block_ * nbf_*nbf_;

	compute_JK();

//	Joe__J->load("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_J.txt");
//	Joe__K->load("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_K.txt");

//	for (int i = 0; i < 24; i++){
//		for (int j = 0; j < 24; j++){
//			jdp[i*24+j] = kdp[i*24+j] = 0.0;
//			jdp[i*24+j] = jjp[i*24+j] - mjp[i*24+j];
//			kdp[i*24+j] = jkp[i*24+j] - mkp[i*24+j];
//		}
//	}

//	jdiff->save("/theoryfs2/ds/obrien/Debug/Psi4/diff_J.txt", false, false, true);
//	kdiff->save("/theoryfs2/ds/obrien/Debug/Psi4/diff_K.txt", false, false, true);

	free(five);
	free(twenty_four);

//	printf("%p\n", (void *) c);
// determine the largest shell size
//	for (size_t shell_iter = 0; shell_iter < Q_shells_ ; shell_iter++ ) {
//		if (auxiliary_->shell(shell_iter).nfunction() > biggest_shell) {
//			biggest_shell = auxiliary_->shell(shell_iter).nfunction();
//		}
//	}
//    fprintf(joe_block, "naux_ is %zu\n", naux_);
//  prepare_Q_blocks();
//std::unique_ptr<double[]> quick_block(new double[naux_* nbf_squared]);
//	double* qb_herald = quick_block.get();
	//for (size_t layer_iter = 0; layer_iter < Q_shells_; layer_iter++){
//		quick_block = std::make_unique<double[]>( auxiliary_->shell(layer_iter).nfunction() * nbf_ * nbf_ );
//		qb_herald = quick_block.get();

//	compute_AO_block_Qpq(0, Q_shells_-1, qb_herald, eri);
		//printf("%f\n", qb_herald[1]);
//	for ( size_t big_P = 0; big_P < naux_; big_P++ ) {
//		for ( size_t little_mu = 0; little_mu < nbf_; little_mu++) {
//			for ( size_t little_nu = 0; little_nu < nbf_; little_nu++ ) {
//				fprintf(joe_block, "%f ", qb_herald[ big_P * nbf_squared + little_mu* nbf_ + little_nu] );
//			}
//			fprintf(joe_block,"\n");
//		}
//		fprintf(joe_block,"\n\n");
//	}
//		quick_block.reset();
	//}
//	printf("It's after this\n");
//	fclose(joe_block);
//    compute_AO_block_Qpq();
}

void DirectDFJK::dgemm_trial(size_t i, size_t j, size_t k) {

	double* A;
	double* B;
	double* C;

	std::unique_ptr<double[]> a(new double[1000]);
	std::unique_ptr<double[]> b(new double[1000]);
	std::unique_ptr<double[]> c(new double[1000]); 

	A = a.get();
	B = b.get();
	C = c.get();

	int a_int = 0;
	int b_int = 0;

/*	//loads A
	for (size_t i_iter = 0; i_iter < i+5; i_iter++) {
		for (size_t k_iter = 0; k_iter < k; k_iter++) {
			A[ i_iter * k + k_iter ] = static_cast<double>(a_int);
			a_int++;
		}
	}

	//loads B
	for (size_t k_iter = 0; k_iter < k; k_iter++) {
		for (size_t j_iter = 0; j_iter < j+5; j_iter++) {
			B[k_iter * j + j_iter ] = static_cast<double>(b_int);
			b_int++;
		}
	}*/

	for (int i = 0; i < 1000; i++) A[i] = B [i] = static_cast<double>(i);

	// C = alpha AB + gamma C
	size_t A_rows = i; // 7/5/19 i = 4
	size_t B_cols = j; // 7/5/19 j = 5
	size_t A_cols_B_rows = k;

	size_t A_cols = k;

	C_DGEMM('N', 'N', 4, 5, 6, 1.0, A, 6, B, 5, 0.0, C, 5);

	FILE* mat_file;
	mat_file = fopen("mat_file.txt", "w");
	for (size_t i_iter = 0; i_iter < i+1; i_iter++) {
		for (size_t j_iter = 0; j_iter < j+1; j_iter++) {
			fprintf(mat_file, "%f ",C[i_iter * (j+1) + j_iter]);
		}
		fprintf(mat_file, "\n");
	}
	fclose(mat_file);

	FILE* c_file = fopen("c_file.txt", "w");
	for (size_t c_iter = 0; c_iter < 25; c_iter++) {
		for (size_t c_it = 0; c_it < 5; c_it++) {
			fprintf(c_file, "%f ", C[(c_iter*5)+c_it]);
		}
		fprintf(c_file, "\n");
	}
	fclose(c_file);

	FILE* a_file = fopen("afile.txt", "w");
	FILE* b_file = fopen("bfile.txt", "w");
	for (int i = 0; i < 1000; i++) {
		fprintf(a_file, "%f  ", A[i]);
		fprintf(b_file, "%f  ", B[i]);
	}
	fclose(a_file);
	fclose(b_file);

}

//syntax copied and adapted from DiskDFJK.
//I decided to forgo copying from MemDFJK because
//that code seems more specialized. 
void DirectDFJK::print_header() const {
	if (print_)	{
		outfile->Printf("==> DirectDFJK: Density-Fitted J/K Matrices <==\n\n");

        outfile->Printf("    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        outfile->Printf("    J tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
		outfile->Printf("    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
		if (do_wK_) outfile->Printf("Don't interpret these results.\nWe don't have w integrals in DirectDFJK\n");
		outfile->Printf("    OpenMP threads:    %11d\n", omp_nthread_);
		outfile->Printf("    Integrals threads: %11d\n", df_ints_num_threads_);
		outfile->Printf("    Memory [MiB]:      %11ld\n", (memory_ * 8L) / (1024L*1024L));
		outfile->Printf("    Algorithm:         %11s\n", "Direct");
        outfile->Printf("    Schwarz Cutoff:    %11.0E\n", cutoff_);
        outfile->Printf("    Fitting Condition: %11.0E\n\n", condition_);

        outfile->Printf("   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level("outfile", print_);
	}
}

// In the long term, I'm not sure what belongs here and what
// belongs in common_init() JSOB

// I think this should call the memory estimate and determine the number of blocks.
// Memory will be allocated locally by the compute_JK() method.

void DirectDFJK::preiterations() {
	printf("Entered preiterations()\n");

	// set_df_ints_num_threads(omp_nthread_);

	//prepare blocks tells us at which indices we start and stop each block.
	// I'm making it a separate function because I've seen evidence that 
	// DFHelper does some resizing, and I want to be able to match that.

	// prepares the coulomb metric. Of course, our calculations depend on the
	//  inverse of the coulomb metric and the (-0.5) power of the metric, so
	//  we will prepare those immediately afterwards
	//prepare_metric();
	
	prepare_metric_power(-1.0);
	prepare_metric_power(-0.5);

	//for what it's worth, we should have the density matrix
	//prepare_D_symm();

}

void sparsity_prep_Qpq() { }

void DirectDFJK::sparsity_prep_pQq(){
	printf("started sparsity prep\n");
	int procs = 1;
#ifdef _OPENMP
	procs = df_ints_num_threads_;
#endif

	// variables to hold data
	double global_max_int = 0.0;

	printf("p_shells_ is %zu and nbf_ is %zu procs is %d\n", p_shells_, nbf_, procs);

	std::unique_ptr<double[]> shel_maxes(new double[p_shells_*p_shells_]);
	std::unique_ptr<double[]> integrals(new double[nbf_*nbf_]);


	for (size_t shel_sq_iter = 0; shel_sq_iter < p_shells_*p_shells_; shel_sq_iter++) shel_maxes[shel_sq_iter] = 0.0;


	for (size_t func_sq_iter = 0; func_sq_iter < nbf_*nbf_; func_sq_iter++) integrals[func_sq_iter] = 0.0;

	// AO machinery
	auto rifactory = std::make_shared<IntegralFactory>(primary_, primary_, primary_, primary_);

	std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);
	std::vector<const double*> buffer(procs);



#pragma omp parallel num_threads(procs)
	{
	printf("parallel execution\n");
	}

#pragma omp parallel num_threads(procs)
    {
    int rank = 0;
#ifdef _OPENMP
    rank = omp_get_thread_num();
#endif
    eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
	buffer[rank] = eri[rank]->buffer();
    }

	// Sparsity member data
	schwarz_shell_mask_pQq_.resize(p_shells_*p_shells_);
	schwarz_func_starts_pQq_.resize(nbf_);
	schwarz_func_map_pQq_.resize(nbf_*nbf_);
	schwarz_func_ints_.resize(nbf_);

	std::vector<size_t> dense_its(nbf_);

	//AO shell metadata
	size_t Shell_M_count;
	size_t Shell_M_start;
	size_t Shell_N_count;
	size_t Shell_N_start;
	double val = 0.0;

#pragma omp parallel for schedule(guided) num_threads(procs) private(Shell_M_count, Shell_M_start, Shell_N_count, Shell_N_start, val) reduction(max:global_max_int)
	for(size_t Shell_M = 0; Shell_M < p_shells_; Shell_M++) {
		int rank = 0;
#ifdef _OPENMP
    	rank = omp_get_thread_num();
#endif
        Shell_M_start = primary_->shell(Shell_M).function_index();
		Shell_M_count = primary_->shell(Shell_M).nfunction();
		for (size_t Shell_N = 0; Shell_N < p_shells_; Shell_N++) {
	        Shell_N_start = primary_->shell(Shell_N).function_index();
			Shell_N_count = primary_->shell(Shell_N).nfunction();
			eri[rank]->compute_shell(Shell_M, Shell_N, Shell_M, Shell_N);
			for (size_t func_m = 0; func_m < Shell_M_count; func_m++){
				for (size_t func_n = 0; func_n < Shell_N_count; func_n++){
					integrals[ (Shell_M_start + func_m) * nbf_ + Shell_N_start + func_n ]
					=
					integrals[ (Shell_N_start + func_n) * nbf_ + Shell_M_start + func_m ]
					=
					val 
					= 
 fabs(buffer[rank][ func_m * Shell_N_count * Shell_M_count * Shell_N_count + 
				    func_n * Shell_M_count * Shell_N_count + 
				    func_m * Shell_N_count + 
				    func_n]);
					if (global_max_int < val) global_max_int = val;
					
					if (shel_maxes[Shell_M * p_shells_ + Shell_N] < val) {
						shel_maxes[Shell_M * p_shells_ + Shell_N]
						=
						shel_maxes[Shell_N * p_shells_ + Shell_M]
						= val;
					}
				}
			}
		}
	}


	double tolerance = cutoff_ * cutoff_ / global_max_int;

	for ( size_t shell_sq_iter = 0; shell_sq_iter < p_shells_; shell_sq_iter++ ) schwarz_shell_mask_pQq_[shell_sq_iter] = ( shel_maxes[shell_sq_iter] < tolerance ? 0 : 1 );

	size_t count;
	for ( size_t func_m = 0; func_m < nbf_; func_m++) {
		count = 0;
		for ( size_t func_n = 0; func_n < nbf_; func_n++ ) {
			if (integrals[ func_m * nbf_ + func_n ] < tolerance) {
				schwarz_func_map_pQq_[func_m * nbf_ + func_n] = 0;
			} else {
				schwarz_func_map_pQq_[func_m * nbf_ + func_n] = count;
				count++;
			}
		}
		schwarz_func_ints_[func_m] = count;
	}

	schwarz_func_starts_pQq_[0] = 0;
	for ( size_t func_it = 1; func_it < nbf_; func_it++ ) {
		schwarz_func_starts_pQq_[func_it] = schwarz_func_starts_pQq_[func_it - 1] + dense_its[func_it - 1] * naux_;
	}
	printf("ended sparsity prep\n");
}

// calculates a block of AO's then conducts the proper contractions to get them
//   where they need to go. I don't think it makes sense to have seperate J and
//   K function calls as J and K are never fully built in a single function.
void DirectDFJK::compute_JK() {
    printf("compute_JK\n");
	if ( pQq_ ) {
    pQp();
	//build_jk_CC_pQq_blocks();
	}  else {

		if (num_blocks_ == 1) {
			printf("Qpq direct\n");
			build_jk_CC_Qpq_direct();
		} else {
			printf("Qpq blocks\n");
    		build_jk_CC_Qpq_blocks();
			printf("Qpq blocks done\n");	
		}

	}

}

void DirectDFJK::postiterations() {
	printf("Entered postiterations\n");
}

// Determines the maximum amount of memory that DirectDFJK will be using
//   at any one time. To be sure, the memory used to construct AO blocks
//   will only be allocated during the compute_JK() method, so it will
//   be available during other procedures. 
// Uses blocking and sparsity information to set total_needs which is 
//   returned by memory_estimate.
void DirectDFJK::our_needs(){
	printf("our needs\n");
	free_memory_ = memory_;
	size_t blocks;

	size_t our_needs = 0;
	double charges_f = 0.0;
	size_t charges_z;
	std::shared_ptr<Molecule> mol_ptr = primary_->molecule();
	
	for (int atom_iter = 0; atom_iter < mol_ptr->natom(); atom_iter++)
		charges_f+=mol_ptr->fZ(atom_iter);

	charges_f = charges_f/2.0;

	charges_z = static_cast<size_t>(charges_f);
	++charges_z;

	//D
	our_needs += nbf_*nbf_;
	free_memory_-=nbf_*nbf_;

	//C_a
	our_needs += nbf_*charges_z;
	free_memory_-=nbf_*charges_z;
	//C_b
	our_needs += nbf_*charges_z;
	free_memory_-=nbf_*charges_z;
	//K
   	our_needs += nbf_*nbf_;	
	free_memory_-=nbf_*nbf_;
   	//J
   	our_needs += nbf_*nbf_;
	free_memory_-=nbf_*nbf_;
	//[ J^{-1} ]
	our_needs += naux_*naux_;
	free_memory_-=naux_*naux_;
	//[ J^{-\frac{1}{2} ]
	our_needs += naux_*naux_;
	free_memory_-=naux_*naux_;
	blocks = (2* (nbf_ * nbf_ *naux_))/free_memory_ +1;
	ABX_block_size_ = ((naux_/blocks)+1)*nbf_*nbf_;

    total_needs_ = our_needs;// + 2*ABX_block_size_;
}
size_t DirectDFJK::memory_estimate(){
	printf("memory_estimate()\n");
	if (pQq_) { 
//		sparsity_prep_pQq();
		prepare_p_blocks();
	}

	if (Qpq_) {
		printf("preparing Q blocks\n");
		if (Qpq_store_sparse_){ sparsity_prep_Qpq(); }
		prepare_Q_blocks();
	}

	our_needs();
	return total_needs_;
}

void DirectDFJK::prepare_metric_power(double power){
	printf("\nmet_power\n");
	SharedMatrix my_met;
	auto coul_met = std::make_shared<FittingMetric>(auxiliary_, true);
	coul_met->form_fitting_metric();
	my_met = coul_met->get_metric();
	my_met->power(power, condition_);
	metric_.push_back(my_met);
	met_powers_.push_back(power);
}

double* DirectDFJK::get_metric_power(double power){
	bool pow_on = false;
	double ret_power;
	size_t met_ind;
	for (size_t met_iter = 0; met_iter < metric_.size(); met_iter++){
		if (std::fabs(power - met_powers_[met_iter]) <condition_ ){
			//ret_power = met_powers_[met_iter];
			met_ind = met_iter;
			pow_on = true;
		}
	}
	if (!pow_on) {
		met_ind = metric_.size();
		SharedMatrix my_met;
		auto coul_met = std::make_shared<FittingMetric>(auxiliary_, true);
		coul_met->form_fitting_metric();
		my_met = coul_met->get_metric();
		my_met->power(power, condition_);
		//prepare_metric_power(power);
		metric_.push_back(my_met);
		met_powers_.push_back(power);
	}
	return metric_[met_ind]->pointer()[0];
}

// called after sparsity_prep_pQq() which is called in the memory estimate
// is called in the memory estimate because the maximum permissable 
// block size is needed for an effective memory estimate, and determining
// that is equivalent to estimating the total amount of memory
// Determines the number of blocks and 
void DirectDFJK::prepare_p_blocks() {
	printf("start prepare_p_blocks()\n");
	size_t current_costs = 0;
	// block max functions
	size_t bmf = 0;
	// atomic charges float
	double charges_f;
	// atomic charges integer
	size_t charges_z;

	size_t block_size = 0;
	size_t shell_funcs;
    size_t biggest_shell = 0;
    
/*	for (size_t func_it = 0; func_it < schwarz_func_ints_.size(); func_it++) {
		sparse_fpf_ += schwarz_func_ints_[func_it];
	}	*/

// We need to figure out the biggest shell in the primary basis set
    // to allocate memory for AO construction
    for (size_t i = 0; i < primary_->nshell(); i++) {
        if (primary_->shell(i).nfunction() > biggest_shell ) {
            biggest_shell = primary_->shell(i).nfunction();
        }
    }
    
	std::shared_ptr<Molecule> mol_ptr = primary_->molecule();

	for (int atom_iter = 0; atom_iter < mol_ptr->natom(); atom_iter++)
		charges_f+=mol_ptr->fZ(atom_iter);
                                                                   
	charges_f = charges_f/2.0;
                                                                   
	charges_z = static_cast<size_t>(charges_f);
	++charges_z;

//
    
// D_ao_ and D_
	current_costs += 3*sizeof(size_t)*nbf_*nbf_;
// Coulomb Metrics
	current_costs += 2*sizeof(double)*naux_*naux_;
	current_costs += 4*sizeof(double)*nbf_*nbf_;
	current_costs += 4*sizeof(double)*nbf_*nbf_;
// adding in costs for U and A
    current_costs += sizeof(double)*biggest_shell*naux_*charges_z;
    current_costs += sizeof(double)*naux_*charges_z;
    
	Shell_starts_.push_back(0);
    
	for (size_t shell_iter = 0; shell_iter < p_shells_; shell_iter++) {
		shell_funcs = primary_->shell(shell_iter).nfunction();
// The 5 is where we're getting our fudge factor we're getting charges_z
//   to stand in for n_occ, and there should only be a factor of 2 included for
//   the number of tensors with dimension occ. The extra 3 is a massive amount of 
//   memory that should cover the overhead from the reast of the calculation.
//   charges_z also includes a fudge factor
		if (  (block_size + shell_funcs + 2) * naux_ * ( 2*charges_z ) >  memory_ - current_costs ){
			Block_funcs_.push_back(block_size);
			block_size = shell_funcs;
			Shell_stops_.push_back(shell_iter - 1);
			Shell_starts_.push_back(shell_iter);
		} else {
			block_size += shell_funcs;
			if (block_size > bmf) {
				 bmf = block_size;
			}
		}
	}
    biggest_block_ = bmf * naux_;
	Shell_stops_.push_back(p_shells_ - 1);

	Block_funcs_.push_back(block_size);

	biggest_block_ = bmf * naux_ * nbf_;
	num_blocks_ = Block_funcs_.size();
    k_disps_.resize(num_blocks_);
    size_t row_disp = 0;
    size_t col_disp = 0;

    for (size_t i = 0; i < num_blocks_; i++) {
        col_disp = 0;
        for (size_t j = 0; j < num_blocks_; j++) {
            k_disps_[i].push_back(row_disp * nbf_ + col_disp);
            col_disp += Block_funcs_[j];
        }
        row_disp+=Block_funcs_[i];
    }

//	schwarz_func_lstarts_pQq_.resize(nbf_);

//	for (size_t i = 0; i < Shell_starts_.size(); i++) {
//		printf("Shell_starts_[%zu] is %zu\n", i, Shell_starts_[i]);
//	}

/*
	size_t subtrahend;
	for (size_t i = 0; i < Shell_starts_.size(); i++) {
		subtrahend = schwarz_func_starts_pQq_[primary_->shell(i).function_index()];
		for (size_t j = 0; j < primary_->shell(i).nfunction(); j++) {
			schwarz_func_lstarts_pQq_[primary_->shell(i).function_index()+j] = schwarz_func_starts_pQq_[primary_->shell(i).function_index()+j] - subtrahend;
		}
	}
 */

	printf("end prepare_p_blocks()\n");
}

//decides where all the blocks will stop and start
// initializes Shell_starts_ and Shell_stops_
// initializes num_blocks_
void DirectDFJK::prepare_Q_blocks() {
//	std::vector<size_t> starts;//Shell_starts_.get();
//	std::vector<size_t> stops;//stops_.get();
	size_t current_costs = 0;
	double charges_f = 0.0;
	size_t charges_z = 0;

	std::shared_ptr<Molecule> mol_ptr = primary_->molecule();

	for (int atom_iter = 0; atom_iter < mol_ptr->natom(); atom_iter++)
		charges_f+=mol_ptr->fZ(atom_iter);
                                                                   
	charges_f = charges_f/2.0;
                                                                   
	charges_z = static_cast<size_t>(charges_f);
	++charges_z;

	current_costs += 3*sizeof(size_t)*nbf_*nbf_;
	current_costs += 2*sizeof(double)*naux_*naux_;
	current_costs += 4*sizeof(double)*nbf_*nbf_;
	current_costs += 4*sizeof(double)*nbf_*nbf_;

	size_t block_max = (free_memory_)/2;
	size_t block_size = 0;
	size_t shell_funcs;

	Shell_starts_.push_back(0);

	biggest_block_ = 0;

	for( size_t shell_iter = 0; shell_iter < Q_shells_; shell_iter++) {
		shell_funcs = auxiliary_->shell(shell_iter).nfunction();
		if ((shell_funcs+block_size+3)*nbf_*(nbf_+2*charges_z)>memory_-current_costs) {
			Block_funcs_.push_back(block_size);
			block_size = shell_funcs;
			Shell_stops_.push_back(shell_iter - 1);//[block_iter] = shell_iter - 1;
			Shell_starts_.push_back(shell_iter);//[block_iter] = shell_iter;
		} else {
			block_size += shell_funcs ;
			if ( block_size *nbf_ * nbf_> biggest_block_ ) { 
				biggest_block_ = block_size * nbf_ * nbf_; 
			}
		}
	}

	Block_funcs_.push_back(block_size);

	Shell_stops_.push_back(Q_shells_-1 );
	
	num_blocks_ = Shell_starts_.size();
}

void DirectDFJK::compute_dense_AO_block_p_pQq(size_t shell, double* ao_block, std::vector<std::shared_ptr<TwoBodyAOInt>> eri){
    int procs = 1;
    
    std::vector<const double*> buffer(procs);
    printf("here from ao\n");

    //int rank = 0;
    double* a  = (double *) malloc(1*sizeof(double));
    
    printf("eri has size %zu\n", eri.size());
    
    buffer[0] = eri[0]->buffer();
    printf("here from ao\n");

    size_t shell_count = primary_->shell(shell).nfunction();
    size_t ShellP_count;
    size_t ShellN_count;
    size_t Shell_ind_0;
    size_t Buff_ind_0;
// the value put in mat_size is nbf*naux
    size_t mat_size = nbf_ * naux_;
    
// In Matt's Code and Rob's code (and, indeed, some of mine), we see 3 loops.
//   however, in this code, we only need two because we only have one shell
//   to compute in the slowest running index, so we're just going to write
//   2 loops to accomodate that.
    printf("shell is %zu p_shells_ is %zu \n", shell, p_shells_);
    
    for (size_t ShellP = 0; ShellP < Q_shells_; ShellP++) {
        ShellP_count = primary_->shell(ShellP).nfunction();
        for (size_t ShellN = 0; ShellN < p_shells_; ShellN++ ) {
            ShellN_count = primary_->shell(ShellN).nfunction();
            Shell_ind_0 = ShellP_count * ShellN_count;
            Buff_ind_0 = ShellN_count * shell_count;
            printf("shell is %zu ShellP is %zu ShellN is %zu\n", shell, ShellP, ShellN);
            eri[0]->compute_shell( ShellP, 0, shell, ShellN);
            int i = 0;
            for (size_t func_m  = 0; func_m < shell_count; func_m++){
                for (size_t func_p = 0; func_p < ShellP_count; func_p++ ) {
                    for (size_t func_n = 0; func_n < ShellN_count; func_n++) {
                        i++;
                        printf("%d\n", i);
                        ao_block[func_m * Shell_ind_0 + func_p * ShellN_count + func_n]
                        a[0] =
                        buffer[0][func_p * Buff_ind_0 + func_m * ShellN_count + func_n];
                    }
                }
            }
        }
    }
    printf("over from ao\n");
}
    
    
void DirectDFJK::compute_sparse_AO_block_p_pQq(size_t start_p, size_t stop_p, double* ao_block, std::vector<std::shared_ptr<TwoBodyAOInt>> eri){
	int procs = 1;
#ifdef _OPENMP
	procs = df_ints_num_threads_;
#endif
    
	std::vector<const double*> buffer(procs);
#pragma omp parallel num_threads(procs)
	{
	int rank = 0;
#ifdef _OPENMP
	rank = omp_get_thread_num();
#endif
	buffer[rank] = eri[rank]->buffer();
	}

// ERI indexing Variables
// Incides are done with romanization of greek letters \Mu -> M \Nu ->N
	size_t ShellP_start;
    size_t ShellP_count;
	size_t ShellM_start;
	size_t ShellM_count;
	size_t ShellN_start;
	size_t ShellN_count;
	size_t Shell_ind_0;

// If we try to save something at &ao_block + ShellP_start *nbv_^2, we
//   will get a segfault in all but the 0th block. disp short for displacemtnt
	size_t ShellP_disp = auxiliary_->shell(start_p).function_index();

	size_t nbf_naux_ = nbf_*naux_;

// We're building A_\mu P \nu
// Outer Loop over prims
#pragma omp parallel for schedule(guided) num_threads(procs) private(ShellP_start, ShellP_count, ShellM_start, ShellM_count, ShellN_start, ShellN_count, Shell_ind_0)
	for ( size_t ShellM = start_p; ShellM <= stop_p; ShellM++ ) {
		int rank = 0;
#ifdef _OPENMP
		rank = omp_get_thread_num();
#endif

		ShellM_start = primary_->shell(ShellM).function_index();
		ShellM_count = primary_->shell(ShellM).nfunction();
// Loop over auxes
		for ( size_t ShellP = start_p; ShellP <= stop_p; ShellP++ ) {
    		ShellP_start = auxiliary_->shell(ShellP).function_index();
    		ShellP_count = auxiliary_->shell(ShellP).nfunction();
// Inner Loop over prims
			for ( size_t ShellN = 0; ShellN < p_shells_; ShellN++ ) {
				if (!schwarz_shell_mask_pQq_[ShellM*p_shells_ + ShellN]){
					continue;
				}
				ShellN_start = primary_->shell(ShellN).function_index();
				ShellN_count = primary_->shell(ShellN).nfunction();
				Shell_ind_0 = ShellM_count*ShellN_count;
// Calculate AO's
				eri[rank]->compute_shell(ShellP, 0, ShellM, ShellN);
// repack AO's
				for( size_t intp = 0; intp < ShellP_count; intp++ ){
					for( size_t intm = 0; intm < ShellM_count; intm++ ){
						for( size_t intn = 0; intn < ShellN_count; intn++ ){
							if (!schwarz_func_map_pQq_[ intm * nbf_ + intn]) {
								continue;
							}
	ao_block[ schwarz_func_starts_pQq_[intm] + schwarz_func_ints_[intm] * intp + intn] 
									=
	//ao_block[ schwarz_func_starts_pQq_[intn] + schwarz_func_ints_[intn] * intp + intm]
    //                        		=
							buffer[rank][ intp * Shell_ind_0  +
                                       intm * ShellN_count +
                                       intn ]; 
						}// inner loop over prim ints
					}// outer loop over prim ints
				}// Loop over aux ints
			}// Inner Loop over prim shells
		}// Outer Loop over prim shells
	}// Loop over aux shells 
}

//basic block AO function to get myself off the ground. 
//this might not be the optimal blocking scheme to calculate AO's
//I'll feel just fine nuking this and testing replacements as the
//project progresses. JSOB

// A note on functionality: it might seem odd that we're passing member data to
//   our own function, but we don't want this function to think that hard.
//   Passing the ao block might also seem odd, but this will be the biggest 
//   single piece of memory used in the scf calculation, so we don't want to
//   reallocate it all the time, and we can't keep it around as member data 
//   because then it will live forever, and that's too much memory to keep
//   around. The memory in &ao_block will live in a unique pointer in the 
//   module that wraps this function, and in this function, we'll just worry 
//   about filling it with what belongs inside.
void DirectDFJK::compute_AO_block_Qpq(size_t start_Q, size_t stop_Q, double* ao_block, std::vector<std::shared_ptr<TwoBodyAOInt>> eri){
	int procs = 1;
#ifdef _OPENMP
	procs = df_ints_num_threads_;
#endif

	std::vector<const double*> buffer(1);

#pragma omp parallel num_threads(procs)
	{
	int rank = 0;
#ifdef _OPENMP
	rank = omp_get_thread_num();
#endif
	buffer[rank] = eri[rank]->buffer();
	}


 // ERI indexing Variables
 // Incides are done with romanization of greek letters \Mu -> M \Nu ->N
  	size_t ShellP_start;
    size_t ShellP_count;
 	size_t ShellM_start;
 	size_t ShellM_count;
 	size_t ShellN_start;
 	size_t ShellN_count;
 	size_t Shell_ind_0; 

// If we try to save something at &ao_block + ShellP_start *nbv_^2, we
//   will get a segfault in all but the 0th block. disp short for displacemtnt
	size_t ShellP_disp = auxiliary_->shell(start_Q).function_index();

	size_t nbf_squared = nbf_*nbf_;


// We're building A^P_{\mu \nu}
// Loop over auxes
	timer_on("DDF AO_CONST");
#pragma omp parallel for schedule(guided) num_threads(procs) private(ShellP_start, ShellP_count, ShellM_start, ShellM_count, ShellN_start, ShellN_count, Shell_ind_0)
	for ( size_t ShellP = start_Q; ShellP <= stop_Q; ShellP++ ) {
		int rank = 0;
#ifdef _OPENMP
		rank = omp_get_thread_num();
#endif
		ShellP_start = auxiliary_->shell(ShellP).function_index();
		ShellP_count = auxiliary_->shell(ShellP).nfunction();
// Outer Loop over prims
		for ( size_t ShellM = 0; ShellM < p_shells_; ShellM++ ) {
			 ShellM_start = primary_->shell(ShellM).function_index();
			 ShellM_count = primary_->shell(ShellM).nfunction();
// Inner Loop over prims
			for ( size_t ShellN = ShellM; ShellN < p_shells_; ShellN++ ) {
				 ShellN_start = primary_->shell(ShellN).function_index();
				 ShellN_count = primary_->shell(ShellN).nfunction();
				 Shell_ind_0 = ShellM_count*ShellN_count;
// Calculate AO's
				eri[rank]->compute_shell(ShellP, 0, ShellM, ShellN);
// repack AO's
				for( size_t intp = 0; intp <  ShellP_count; intp++ ){
					for( size_t intm = 0; intm < ShellM_count; intm++ ){
						for( size_t intn = 0; intn < ShellN_count; intn++ ){
							ao_block[ (ShellP_start - ShellP_disp + intp) * nbf_squared + 
									  (ShellM_start + intm)               * nbf_ +
									   ShellN_start + intn ] =
							ao_block[ (ShellP_start - ShellP_disp + intp) * nbf_squared + 
									  (ShellN_start + intn)               * nbf_ +
									   ShellM_start + intm ] =

							buffer[rank][ intp * Shell_ind_0 + 
                                      intm * ShellN_count +
                                      intn ]; 


						}// inner loop over prim ints
					}// outer loop over prim ints
				}// Loop over aux ints
			}// Inner Loop over prim shells
		}// Outer Loop over prim shells
	}// Loop over aux shells 
	timer_off("DDF AO_CONST");
}

// The parameter ind says at which index of the D_ao_ vector we're going to
//   store our new D_ao_.
void DirectDFJK::compute_D_ao(size_t ind) {
	if (ind +1 > D_ao_.size()) {
		while ( D_ao_.size() < ind ) {
			SharedMatrix k = std::make_shared<Matrix>();
			D_ao_.push_back(k);
		}
	} else {
		D_ao_.erase(D_ao_.begin() + ind);
	}
	std::stringstream D_name;
	D_name << "D " << ind << " (AO)";
	D_ao_.insert(D_ao_.begin() + ind, std::make_shared<Matrix>(D_name.str(), C_left_ao_[ind]->rowspi(), C_right_ao_[ind]->rowspi() ));

	double** Dp = D_ao_[ind]->pointer(0);
	double** Clp = C_left_ao_[ind]->pointer(0);
	double** Crp = C_right_ao_[ind]->pointer(0);

	int cl_rows = C_left_ao_[ind]->nrow();
	int cl_cols = C_left_ao_[ind]->ncol();
	
	int cr_rows = C_right_ao_[ind]->nrow();
	int cr_cols = C_right_ao_[ind]->ncol();

	C_DGEMM('N', 'T', cl_rows, cr_rows, cl_cols, 1.0, Clp[0], cl_cols, Crp[0], cr_cols, 0.0, Dp[0], cr_rows);
//  C_DGEMM('N', 'T', C_left_ao_[ind]->rowspi(0)[0], C_right_ao_[ind]->rowspi(0)[0], C_left_ao_[ind]->colspi(0)[0], 1.0, Clp[0], C_left_ao_[ind]->colspi(0)[0], Crp[0], C_right_ao_[ind]->colspi(0)[0], 0.0, Dp[0], C_left_ao_[ind] );
	
}

//start and :
void DirectDFJK::V_gets_AD(size_t stop, double* v, double* a, double* d) {
	size_t nbf_squared = nbf_ * nbf_;

	for (size_t v_iter = 0; v_iter < stop; v_iter++) {
		double dub_sum = 0.0;
		for (size_t mu_iter = 0; mu_iter < nbf_; mu_iter++) {
			for (size_t nu_iter = 0; nu_iter < nbf_; nu_iter++) {
				dub_sum += a[ v_iter * nbf_squared + mu_iter * nbf_ + nu_iter ] * d[mu_iter * nbf_ + nu_iter];
			}
		}
		v[v_iter] = dub_sum;
	}
}

void DirectDFJK::Accumulate_J(size_t stop, double* j, double* a, double* phi){
	size_t nbf_squared = nbf_ * nbf_;
	for (size_t slice_iter = 0; slice_iter < stop; slice_iter++) {
		for (size_t mu_iter = 0; mu_iter < nbf_; mu_iter++) {
			for (size_t nu_iter = 0; nu_iter < nbf_; nu_iter++) {
				j[ (mu_iter * nbf_) + nu_iter] += phi[slice_iter]*a[slice_iter * nbf_squared + mu_iter * nbf_ + nu_iter];
			}
		}
	}
}


void DirectDFJK::U_gets_AC(size_t stop, size_t c_cols, double* u, double* a, double* c){
	size_t nbf_squared = nbf_ * nbf_;
	size_t u_slice_size = c_cols * nbf_;

	for ( size_t aux_slice_iter = 0; aux_slice_iter < stop; aux_slice_iter++ ) {
		C_DGEMM('N', 'N', nbf_, c_cols, nbf_, 1.0, &a[aux_slice_iter * nbf_squared], nbf_, &c[0], c_cols, 0.0, &u[ aux_slice_iter * u_slice_size ], c_cols);
	}

/*
	double nu_a_sum;
	for (size_t aux_slice_iter = 0; aux_slice_iter < stop; aux_slice_iter++){
		for (size_t nu_iter = 0; nu_iter < nbf_; nu_iter++ ){
			for (size_t a_iter = 0; a_iter < c_cols; a_iter++){
				nu_a_sum = 0;
				for (size_t sigma_iter = 0; sigma_iter < nbf_; sigma_iter++){
					nu_a_sum += a[aux_slice_iter*nbf_squared + nu_iter*nbf_ + sigma_iter] * c[sigma_iter*c_cols + a_iter];
				}
				u[aux_slice_iter*u_slice_size + nu_iter*c_cols + a_iter] = nu_a_sum;
			}
		}
	}
*/
}

// We will assume that X comes displaced for the appropriate block
void DirectDFJK::X_accumulates_JU(size_t x_stop, size_t u_stop, size_t slice_size, size_t x_off, size_t u_off, size_t c_cols, double* x, double* j, double* u) {
	double* our_j;
	double met_factor;
	for (size_t x_slice = 0; x_slice < x_stop; x_slice++){
		for (size_t u_slice = 0; u_slice < u_stop; u_slice++) {
			met_factor = j[(u_slice+u_off)*naux_ + x_off + x_slice];
			for (size_t nu_iter = 0; nu_iter < nbf_; nu_iter++) {
				for (size_t a_iter = 0; a_iter < c_cols; a_iter++) {
					x[x_slice*slice_size + nu_iter * c_cols + a_iter ] += met_factor * u[u_slice*slice_size + nu_iter *c_cols + a_iter] ;
				}
			}
		}
	}
}

void DirectDFJK::Accumulate_K_c_is_c( size_t stop, size_t slice_size, size_t c_cols, double* k, double* x){

	for ( size_t slice_iter = 0; slice_iter < stop; slice_iter++ ) {
		C_DGEMM('N', 'T', nbf_, nbf_, c_cols, 1.0, &x[slice_size * slice_iter] , c_cols, &x[slice_size * slice_iter], c_cols, 1.0, k, nbf_);
	}

/*
	double mu_nu_sum;
	for (size_t slice_iter = 0; slice_iter < stop; slice_iter++) {
		for (size_t mu_iter = 0; mu_iter < nbf_; mu_iter++) {
			for (size_t nu_iter = 0; nu_iter < nbf_; nu_iter++) {
				mu_nu_sum = 0.0;
				for (size_t a_iter = 0; a_iter < c_cols; a_iter++){
					mu_nu_sum += x[slice_iter*slice_size + mu_iter*c_cols + a_iter] * x[slice_iter*slice_size + nu_iter*c_cols + a_iter];
				}
				k[mu_iter*nbf_ + nu_iter ] += mu_nu_sum;
			}
		}
	}
*/

}

void DirectDFJK::build_jk_CC_Qpq_direct() {
    int procs = 1;
#ifdef _OPENMP
	procs = df_ints_num_threads_;
#endif

	std::shared_ptr<BasisSet> zero_ = BasisSet::zero_ao_basis_set();
	auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero_, primary_, primary_);
	std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);


#pragma omp parallel num_threads(procs)
	{
	int rank = 0;
#ifdef _OPENMP
	rank = omp_get_thread_num();
#endif
	eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
	}

	x_size_ = C_left_ao_[0]->ncol() * nbf_ * naux_;

	double* c = C_left_ao_[0]->pointer()[0];
	double* d = D_ao_[0]->pointer()[0];
	double* met_m_1_0 = get_metric_power(-1.0);
	double* met_m_0_5 = get_metric_power(-0.5);

	printf("%d %zu %zu \n", C_left_ao_[0]->ncol() , nbf_ , naux_);


	std::unique_ptr<double[]> A(new double[nbf_*nbf_*naux_]);
	std::unique_ptr<double[]> U(new double[x_size_]);
	std::unique_ptr<double[]> X(new double[x_size_]);
	std::unique_ptr<double[]> V(new double[naux_]);
	std::unique_ptr<double[]> PHI(new double[naux_]);

	double* a = A.get();
	double* u = U.get();
	double* x = X.get();
	double* v = V.get();
	double* phi = PHI.get();

	double* j = J_ao_[0]->pointer()[0];
	double* k = K_ao_[0]->pointer()[0];

	x_slice_.push_back(nbf_*C_left_ao_[0]->ncol());

	size_t nbf_squared = nbf_* nbf_;

	for(size_t x_iter = 0; x_iter < x_size_; x_iter++){
		x[x_iter] = 0.0;
	}

	//We'll be accumulating into both these vectors. Make them zero.
	for(size_t phi_v_iter = 0; phi_v_iter < naux_; phi_v_iter++){
		phi[phi_v_iter] = 0.0;
		v[phi_v_iter] = 0.0;
	}

	compute_AO_block_Qpq( Shell_starts_[0], Shell_stops_[0], a, eri);

//	V_gets_AD( naux_, v + in_block_off, a, d);
	C_DGEMV( 'N', (int) naux_, (int) nbf_squared, 1.0,         a, (int) nbf_squared,   d, 1, 0.0,   v, 1);

	// PHI /gets \cmpq V
	C_DGEMV( 'T', (int) naux_, (int)       naux_, 1.0, met_m_1_0, (int)       naux_,   v, 1, 0.0, phi, 1);

//	Accumulate_J( Block_funcs_[0] , j, a, phi);
// Accumulate_J
	C_DGEMV( 'T', (int) naux_, (int) nbf_squared, 1.0,         a, (int) nbf_squared, phi, 1, 0.0,   j, 1);

//	U_gets_AC(Block_funcs_[0], static_cast<size_t>(C_left_ao_[0]->ncol()), u, a, c);
// Form U by flattening out U over the auxiliary basis.
	C_DGEMM( 'N', 'N', (int) (nbf_*naux_), C_left_ao_[0]->ncol(), (int) nbf_, 1.0, a, (int) nbf_, c, C_left_ao_[0]->ncol(), 0.0, u, C_left_ao_[0]->ncol());

//	X_accumulates_JU(Block_funcs_[0], Block_funcs_[0], x_slice_[0], ou_block_off, in_block_off, static_cast<size_t>(C_left_ao_[0]->ncol()), x, met_m_0_5, u );
// Form x by flattening out U along the auxiliary basis and transposing
// to contract over the auxiliary basis set.
	C_DGEMM( 'T', 'T', ((int) nbf_) * C_left_ao_[0]->ncol(), (int) naux_, (int) naux_, 1.0, u, ((int) nbf_) * C_left_ao_[0]->ncol(), met_m_0_5, (int) naux_, 0.0, x, (int) naux_);

//	Accumulate_K_c_is_c(Block_funcs_[0], x_slice_[0], static_cast<size_t>(C_left_ao_[0]->ncol()), k, x);
	timer_on("Big DGEMM");
	C_DGEMM( 'N', 'T', (int) nbf_, (int) nbf_, ((int) naux_) * C_left_ao_[0]->ncol(), 1.0, x, ((int) naux_) * C_left_ao_[0]->ncol(), x, ((int) naux_) * C_left_ao_[0]->ncol()  , 0.0, k, (int) nbf_);
	timer_off("Big DGEMM");

//	J_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_J.txt", false, false, true);
//	K_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_K.txt", false, false, true);

//	C_left_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/C_directdfjk.txt", false, false, true);

	x_slice_.erase(x_slice_.begin(), x_slice_.end() );
}


void DirectDFJK::build_jk_CC_Qpq_blocks() {
    int procs = 1;
#ifdef _OPENMP
	procs = df_ints_num_threads_;
#endif

	std::shared_ptr<BasisSet> zero_ = BasisSet::zero_ao_basis_set();
	auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero_, primary_, primary_);
	std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);

#pragma omp parallel num_threads(procs)
	{
	int rank = 0;
#ifdef _OPENMP
	rank = omp_get_thread_num();
#endif
	eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
	}

	x_size_ = (biggest_block_/nbf_) * C_left_ao_[0]->ncol();

	double* c = C_left_ao_[0]->pointer()[0];
	double* d =  D_ao_[0]->pointer()[0];
	double* met_m_1_0 = get_metric_power(-1.0);
	double* met_m_0_5 = get_metric_power(-0.5);

	std::unique_ptr<double[]> HCMT(new double[naux_*naux_]);

	double* met_m_0_5_T = HCMT.get();

	for (size_t ni = 0; ni < naux_; ni++){
		for (size_t nj = 0; nj < naux_; nj++){
			met_m_0_5_T[ni * naux_ + nj] = met_m_0_5[nj*naux_ + ni];
		}
	}

	printf("nbf is %zu, naux_ is %zu \n", nbf_, naux_ );

	for (size_t i = 0; i < Shell_starts_.size(); i++) {
		printf("Shell_starts_[%zu] is %zu \n", i, Shell_starts_[i]);
		printf("Shell_stops_[%zu] is %zu \n", i, Shell_stops_[i]);
	}

	printf("Shell_starts.size() is %zu Shell_stops_.size() is %zu \n", Shell_starts_.size(), Shell_stops_.size());

	std::unique_ptr<double[]> A(new double[biggest_block_]);
	std::unique_ptr<double[]> U(new double[x_size_]);
	std::unique_ptr<double[]> X(new double[x_size_]);
	std::unique_ptr<double[]> V(new double[naux_]);
	std::unique_ptr<double[]> PHI(new double[naux_]);

	double* a = A.get();
	double* u = U.get();
	double* x = X.get();
	double* v = V.get();
	double* phi = PHI.get();

	double* j = J_ao_[0]->pointer(0)[0];
	double* k = K_ao_[0]->pointer(0)[0];

	for (size_t jk_iter = 0; jk_iter < nbf_*nbf_; jk_iter++) j[0] = k[0] = 0.0;

	x_slice_.push_back(nbf_*C_left_ao_[0]->ncol());

	//Inner Block offset. Used for pointer math.
	size_t in_block_off  {0};
	size_t ou_block_off  {0};

	size_t nbf_squared = nbf_*nbf_;

	for(size_t phi_v_iter = 0; phi_v_iter < naux_; phi_v_iter++){
		phi[phi_v_iter] = 0.0;
		v[phi_v_iter] = 0.0;
	}

	for (size_t block_iter_ou = 0; block_iter_ou < num_blocks_; block_iter_ou++) {
// there should be a blas call for this
		for(size_t x_iter = 0; x_iter < x_size_; x_iter++){
			x[x_iter] = 0.0;
		}
		in_block_off = 0;
		if ( block_iter_ou == 1 ) {
			C_DGEMV('N', naux_, naux_, 1.0, met_m_1_0, naux_, v, 1, 0.0, phi, 1);
		}
		for (size_t block_iter_in = 0; block_iter_in < num_blocks_; block_iter_in++) {
			printf("about to compute aos\n");
			compute_AO_block_Qpq( Shell_starts_[block_iter_in], Shell_stops_[block_iter_in], a, eri);
			printf("computed aos\n");
			if ( block_iter_ou == 0 ) {
//				V_gets_AD( Block_funcs_[block_iter_in], v + in_block_off, a, d);// + v_phi_add);
				C_DGEMV( 'N', Block_funcs_[block_iter_in], nbf_squared, 1.0, a, nbf_squared,   d, 1, 0.0, v + in_block_off, 1);
			}
			if ( block_iter_ou == 1 ) {
//				Accumulate_J( Block_funcs_[block_iter_in] , j, a, phi + in_block_off);
				C_DGEMV( 'T', Block_funcs_[block_iter_in], nbf_squared, 1.0, a, nbf_squared, phi + in_block_off, 1, 1.0, j, 1);
			}
//			U_gets_AC(Block_funcs_[block_iter_in], static_cast<size_t>(C_left_ao_[0]->ncol()), u, a, c);
			C_DGEMM( 'N', 'N', Block_funcs_[block_iter_in] * nbf_, C_left_ao_[0]->ncol(), nbf_, 1.0, a, nbf_, c, C_left_ao_[0]->ncol(), 0.0, u, C_left_ao_[0]->ncol());
//			X_accumulates_JU(Block_funcs_[block_iter_ou], Block_funcs_[block_iter_in], x_slice_[0], ou_block_off, in_block_off, static_cast<size_t>(C_left_ao_[0]->ncol()), x, met_m_0_5, u );
			C_DGEMM( 'T', 'N', nbf_* C_left_ao_[0]->ncol(), Block_funcs_[block_iter_ou], Block_funcs_[block_iter_in], 1.0, u, nbf_* C_left_ao_[0]->ncol(), met_m_0_5 + (in_block_off * naux_) + ou_block_off, naux_, 1.0, x, Block_funcs_[block_iter_ou]);
			in_block_off += Block_funcs_[block_iter_in];
		}
		C_DGEMM( 'N', 'T', nbf_, nbf_, Block_funcs_[block_iter_ou] * C_left_ao_[0]->ncol(), 1.0, x, Block_funcs_[block_iter_ou] * C_left_ao_[0]->ncol(), x, Block_funcs_[block_iter_ou] * C_left_ao_[0]->ncol(), 1.0, k, nbf_ );
//		Accumulate_K_c_is_c(Block_funcs_[block_iter_ou], x_slice_[0], static_cast<size_t>(C_left_ao_[0]->ncol()), k, x);
		ou_block_off += Block_funcs_[block_iter_ou];
//		x_block = x_block + Block_funcs[block_iter_ou] * x_slice_[0];
	}
//J_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_J.txt", false, false, true);
//K_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_K.txt", false, false, true);
	x_slice_.erase( x_slice_.begin(), x_slice_.end() ); 
 }

void DirectDFJK::pQp(){
    int procs = 1;
    int rank = 0;
    
//dirty
    size_t biggest_shell = 0;
    
    for (int i = 0; i < primary_->nshell(); i++){
        if (primary_->shell(i).nfunction() > biggest_shell) {
            biggest_shell = primary_->shell(i).nfunction();
        }
    }
//end dirty

    std::shared_ptr<BasisSet> zero_ = BasisSet::zero_ao_basis_set();
    auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero_, primary_, primary_);
    std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);
    eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
    
    printf("eri has size %zu\n", eri.size());
    printf("auxiliary_ is %s \n", auxiliary_->name().c_str());
    printf("biggest_shell is %zu \n", biggest_shell);
    printf("auxiliary_ has %d shells\n", auxiliary_->nshell());
    int nocc = C_left_ao_[0]->ncol();
    
    double* k = K_ao_[0]->pointer(0)[0];
    
    std::unique_ptr<double[]> A(new double[biggest_shell*naux_*nbf_]);
    std::unique_ptr<double[]> U(new double[naux_*nbf_]);
    std::unique_ptr<double[]> XN(new double[biggest_block_*nocc]);
    std::unique_ptr<double[]> XO(new double[biggest_block_*nocc]);
    std::unique_ptr<double[]> V(new double[naux_]);

    double* a = A.get();
    double* u = U.get();
    double* xn = XN.get();
    double* xo = XO.get();
    double* v = V.get();
    
    char first_char = ( num_blocks_ == 0 ? 'V' : 'B');
    
    X_Block(first_char, true, 0, a, xo, u, v, eri);
    C_DGEMM('N', 'T', Block_funcs_[0], Block_funcs_[0], naux_*nocc, 1.0, xo, naux_*nocc, xo, naux_*nocc, 0.0, k, nbf_);
    
    for (size_t xo_iter = 0; xo_iter < num_blocks_; xo_iter++){
        for (size_t xn_iter = 1; xn_iter < num_blocks_ - xo_iter; xn_iter++){
            if ( xo_iter == 0 && xn_iter != num_blocks_ - 1) {
                X_Block('V', true, xn_iter, a, xn, u, v, eri);
                C_DGEMM('N', 'T', Block_funcs_[xn_iter], Block_funcs_[xn_iter], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 0.0, k + k_disps_[xn_iter][xn_iter], nbf_ );
                C_DGEMM('N', 'T', Block_funcs_[xo_iter], Block_funcs_[xn_iter], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 0.0, k + k_disps_[xo_iter][xn_iter], nbf_ );
            } else if (xo_iter == 0 && xn_iter == num_blocks_ - 1) {
                X_Block('B', true, xn_iter, a, xn, u, v, eri);
                C_DGEMM('N', 'T', Block_funcs_[xn_iter], Block_funcs_[xn_iter], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 0.0, k + k_disps_[xn_iter][xn_iter], nbf_);
            } else if ( xo_iter != 0 ) {
                X_Block('P', true, xo_iter, a, xn, u, v, eri);
                C_DGEMM('N', 'T', Block_funcs_[xn_iter], Block_funcs_[xo_iter], naux_*nocc, 1.0, xn, naux_*nocc, xo, naux_*nocc, 0.0, k + k_disps_[xn_iter][xo_iter], nbf_);
            }
        }
        xo = xn;
    }
    
    printf("ping 1357\n");
    
    if (first_char != 'B') {
        X_Block('P', false, num_blocks_ - 1, a, nullptr, u, v, eri);
    }
    
    
}

void DirectDFJK::build_jk_CC_pQq_blocks(){
	printf("in build_jk\n");
	int procs = 1;
#ifdef _OPENMP
	procs = df_ints_num_threads_;
#endif

	std::shared_ptr<BasisSet> zero_ = BasisSet::zero_ao_basis_set();
	auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero_, primary_, primary_);
	std::vector<std::shared_ptr<TwoBodyAOInt>> eri(procs);

#pragma omp parallel num_threads(procs)
	{
	int rank = 0;
#ifdef _OPENMP
	rank = omp_get_thread_num();
#endif
	eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
	}

	double* c = C_left_ao_[0]->pointer()[0];

	int nocc = C_left_ao_[0]->ncol();

	double* d = D_ao_[0]->pointer()[0];

	double* met_m_1_0 = get_metric_power(-1.0);
	double* met_m_0_5 = get_metric_power(-0.5);

	double* j = J_ao_[0]->pointer(0)[0];
	double* k = K_ao_[0]->pointer(0)[0];

// some would say too much C not enough C++. I say I like my 
//    pointers to be special. Also, I'm definitely going to 
//    test if malloc makes it go faster.
	double* pruned_j;
	std::unique_ptr<double[]> P_J(new double[nbf_]);
	pruned_j = P_J.get();

	x_size_ = ( biggest_block_ / sparse_fpf_ ) * C_left_ao_[0]->ncol();

	double* pruned_c;
	double* pruned_d;
 	double* a;
	double* xl;
	double* xr;
	double* v;
	double* phi;
	double* ub;


	std::unique_ptr<double[]> P_C(new double[nbf_ * C_left_ao_[0]->ncol()]);
	std::unique_ptr<double[]> P_D(new double[nbf_]);
	std::unique_ptr<double[]> A(new double[biggest_block_]);
	std::unique_ptr<double[]> XL(new double[x_size_]);
	std::unique_ptr<double[]> XR(new double[x_size_]);
	std::unique_ptr<double[]> V(new double[naux_]);
    std::unique_ptr<double[]> PHI(new double[naux_]);
	std::unique_ptr<double[]> UB(new double[nocc*naux_]);

	pruned_c = P_C.get();
	pruned_d = P_D.get();
	a = A.get(); size_t a_cols;
	xl = XL.get(); size_t xl_ind;
	xr = XR.get(); size_t xr_ind;
	v = V.get();
	phi = PHI.get();
	ub = UB.get();

	for (size_t zeroer = 0; zeroer < naux_; zeroer++) {
		v[zeroer] = 0.0;
		phi[zeroer] = 0.0;
	}

printf("line number \n");
// Calculate eri's
	compute_sparse_AO_block_p_pQq( Shell_starts_[0], Shell_stops_[0], a, eri);
	xl_ind = 0;
	for ( size_t nu = 0; nu < Block_funcs_[0] ; nu++ ) {
		prune_pQq(nu, static_cast<size_t>(C_left_ao_[0]->ncol()), pruned_c, c);
// U <- AC one matrix of U and A at a time
printf("line number 1336\n");
		C_DGEMM('N', 'N', naux_, nocc, schwarz_func_ints_[nu], 1.0, a + schwarz_func_lstarts_pQq_[nu], schwarz_func_ints_[nu], pruned_c, nocc, 0.0, ub, nocc );
// X <- \cmpq U
printf("line number 1339\n");
		C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, ub, nocc, 0.0, xl + nu*naux_*nocc, nocc);
	}
// K <- XX
printf("line number 1313\n");
	C_DGEMM('N', 'T', Block_funcs_[0], Block_funcs_[0], naux_*nocc, 1.0, xl, naux_*nocc, xl, naux_*nocc, 0.0, k, nbf_);


// //
	size_t gemvcols = 0;
	for ( size_t lambda = 0; lambda < Block_funcs_[0]; lambda++ ) {
		gemvcols += schwarz_func_ints_[lambda];
	}
	for ( size_t lambda = 0; lambda < Block_funcs_[0]; lambda++ ) {
		prune_pQq(lambda, 1, pruned_d, d + lambda*nbf_ );
// V <- AD
		C_DGEMV('N', naux_, schwarz_func_ints_[lambda], 0.0, a + schwarz_func_lstarts_pQq_[lambda], gemvcols, pruned_d, 1, 1.0, v, 1);
	}



	//Major loops to calculate J & K
	for ( size_t l_iter = 0; l_iter < num_blocks_ - 1; l_iter++ ) {
		
		for ( size_t r_iter = 1; r_iter < num_blocks_ - l_iter; r_iter++) {
			compute_sparse_AO_block_p_pQq( Shell_starts_[r_iter], Shell_stops_[r_iter], a, eri);
// compute V for coulomb matrix construction
			gemvcols = 0;
			for ( size_t lambda = primary_->shell(r_iter).function_index(); lambda < primary_->shell(r_iter).function_index() + Block_funcs_[r_iter]; lambda++ ) {
				gemvcols += schwarz_func_ints_[lambda];
			}
			for ( size_t lambda = primary_->shell(r_iter).function_index(); lambda < primary_->shell(r_iter).function_index() + Block_funcs_[r_iter]; lambda++) {
				prune_pQq(nbf_, 1, pruned_d, d + lambda*nbf_);
				C_DGEMV('T', naux_, schwarz_func_ints_[lambda], 0.0, a + schwarz_func_lstarts_pQq_[lambda], gemvcols, pruned_d, 1, 1.0, v, 1);
			}


			xr_ind = r_iter;
// contract A and C into U , U and cmpq into X
			for ( size_t nu = 0; nu < Block_funcs_[r_iter] ; nu++ ) {
				prune_pQq(nu, static_cast<size_t>(C_left_ao_[0]->ncol()), pruned_c, c);
printf("line number 1350\n");
				C_DGEMM('N', 'N', naux_, nocc, schwarz_func_ints_[nu], 1.0, a + schwarz_func_lstarts_pQq_[nu], schwarz_func_ints_[nu], pruned_c, nocc, 0.0, ub, nocc );
printf("line number 1352\n");
				C_DGEMM('N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, ub, nocc, 0.0, xr + nu*nocc*naux_, nocc);
			} 
// Actually fill elements of k. 
// The i == 0 case fillsin the 0th row and the diagonal blocks of k as:
/*  x x x x x
 *    x
 *      x
 *        x
 *          x
*/
			if ( r_iter == 0) {
printf("line number 1364\n");
				C_DGEMM('N', 'T', Block_funcs_[xl_ind], Block_funcs_[xr_ind], nocc*naux_, 0.0, xl, nocc*naux_, xr, nocc*naux_, 0.0, k + Block_funcs_[xl_ind]*xr_ind,nbf_);
				size_t kdisp = 0;
				for (size_t bl_iter = 0; bl_iter < xr_ind; bl_iter++) {
					kdisp += Block_funcs_[bl_iter] * nbf_;
				}
				kdisp += primary_->shell(Shell_starts_[r_iter]).function_index();
printf("line number 1371\n");
				C_DGEMM('N', 'T', Block_funcs_[xr_ind], Block_funcs_[xr_ind], naux_*nocc, 1.0, xr, naux_*nocc, xr, naux_*nocc, 0.0, k + kdisp, nbf_);
// DAXBY phi = cmpq v;  
				C_DGEMV( 'T', naux_, naux_, 0.0, met_m_1_0, naux_, v, 1, 0.0, phi, 1);
// contract j = A_ mu block nb _ P _ nu sparse mu , phi
				gemvcols = 0;
	        	for ( size_t mu = primary_->shell(xr_ind).function_index(); mu < primary_->shell(xr_ind).function_index() + Block_funcs_[xr_ind]; mu++ ) {
            		gemvcols += schwarz_func_ints_[mu];
            	}
				for (size_t mu = primary_->shell(Shell_starts_[xr_ind]).function_index(); mu < primary_->shell(Shell_starts_[xr_ind]).function_index() + Block_funcs_[xr_ind]; mu++){
					C_DGEMV( 'T', naux_, schwarz_func_ints_[mu], 0.0, a + schwarz_func_lstarts_pQq_[mu], gemvcols, phi, 1, 0.0, pruned_j, 1);
					unprune_j_pQq(mu, pruned_j, j);
				}
			
			} else{
// The else case sees us filling in the (nb - 1)th column of the matrix
//   e.g the blocks when i = 1 and j = 1
/*  x x x x x
 *    x     x
 *      x  
 *        x
 *          x
*/
			size_t kdisp = 0;
			for (size_t bl_iter = 0; bl_iter < xr_ind; bl_iter++) {
				kdisp += Block_funcs_[bl_iter] * nbf_;
			}
			kdisp += primary_->shell(Shell_starts_[xl_ind]).function_index();
printf("line number 1399\n");
			C_DGEMM('N', 'T', Block_funcs_[xl_ind], Block_funcs_[xr_ind], nocc*naux_, 0.0, xl, nocc * naux_, xr, naux_ * nocc, 0.0, k, nbf_);
			}
		}
		xl = xr;
		xl_ind = xr_ind;
	}

	if (num_blocks_) {
		compute_sparse_AO_block_p_pQq( Shell_starts_[0], Shell_stops_[0], a, eri);
		xr_ind  = 0;
	}

	gemvcols = 0;
    for ( size_t mu = primary_->shell(xr_ind).function_index(); mu < primary_->shell(xr_ind).function_index() + Block_funcs_[xr_ind]; mu++ ) {
    	gemvcols += schwarz_func_ints_[mu];
    }
    for (size_t mu = primary_->shell(Shell_starts_[xr_ind]).function_index(); mu < primary_->shell(Shell_starts_[xr_ind]).function_index() + Block_funcs_[xr_ind]; mu++){
    	C_DGEMV( 'T', naux_, schwarz_func_ints_[mu], 0.0, a + schwarz_func_lstarts_pQq_[mu], gemvcols, phi, 1, 0.0, pruned_j, 1);
    	unprune_j_pQq(mu, pruned_j, j);
    }

	for ( size_t l_iter = 0; l_iter < nbf_; l_iter++ ){
		for (size_t r_iter = l_iter ; r_iter < nbf_; r_iter++ ) {
			k[ r_iter * naux_ + l_iter ] = k[l_iter*naux_ + r_iter];
		}
	}

}

void DirectDFJK::prune_pQq( size_t bf, size_t nocc, double* pruned_c, double* raw_c ) {
	for (size_t z_iter = 0; z_iter < nbf_ * nocc; z_iter++) { pruned_c[z_iter] = 0.0; }

	for (size_t row_it = 0; row_it < nbf_; row_it++) {
		if ( schwarz_func_map_pQq_[ bf * nbf_ + row_it] ) {
			for ( size_t col_it = 0; col_it < nocc; col_it++) {
				pruned_c[schwarz_func_map_pQq_[row_it]*nocc + col_it] = raw_c[row_it * nocc + col_it];
			}
		}
	}

}

// Function that produces tensors for the final contraction
//    of the exchange matrix build. The issue is that it may or may not
//    have to construct one of various terms for the coulomb matrix.
//    We handle this with a switch.
// coul_work \in { 'V', 'P', 'B' }
// 'V' means we compute a vector to be contracted against the coulomb Metric.
void DirectDFJK::X_Block( char coul_work, bool compute_k, size_t block, double* ao_block, double* x, double* u, double* coulomb_vector, std::vector<std::shared_ptr<TwoBodyAOInt>> eri){

    printf("pingX coul_work is %c\n", coul_work);
    
    size_t nocc = C_left_ao_[0]->ncol();
    double* c = C_left_ao_[0]->pointer()[0];
    double* met_m_0_5 = get_metric_power(-0.5);
    double* j = J_ao_[0]->pointer()[0];
    double* d = D_ao_[0]->pointer()[0];

    
    switch (coul_work) {
        case 'V':
            for (size_t shell_iter = Shell_starts_[shell_iter]; shell_iter <= Shell_stops_[shell_iter]; shell_iter++){
                //compute ao blocks
                compute_dense_AO_block_p_pQq(shell_iter, ao_block , eri);
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
// Form V for Coulomb Matrix construction
                    C_DGEMV( 'N', (int) naux_, (int) nbf_, 1.0, ao_block + func_it*naux_*nbf_, nbf_, d + primary_->shell(shell_iter).function_index()*nbf_ + func_it * nbf_, 1, 1.0, coulomb_vector, 1 );
                }
                if (compute_k) {
// Form U for Exchange Matrix construction
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, nbf_, 1.0, ao_block, nbf_, c, nocc, 0.0, u, nocc);
// Contract this u into the corresponding portion of x
                    C_DGEMM('N', 'N', naux_, primary_->shell(shell_iter).nfunction()*nocc, naux_, 1.0, met_m_0_5, naux_, u, naux_ * primary_->shell(shell_iter).nfunction(), 0.0, x + naux_ * nbf_ * (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index()), naux_ * primary_->shell(shell_iter).nfunction() );
                }
            }
            break;
        case 'P':
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
                compute_dense_AO_block_p_pQq(shell_iter, ao_block, eri);
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
                    C_DGEMV('T', (int) naux_, (int) nbf_, 1.0, ao_block + func_it * naux_ * nbf_ , nbf_, coulomb_vector, 1, 1.0,  j + naux_*nbf_ * (primary_->shell(shell_iter).function_index() + func_it ), 1 );
                }
                if (compute_k) {
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, nbf_, 1.0, ao_block, nbf_, c, nocc, 0.0, u, nocc);
                    C_DGEMM('N', 'N', naux_, primary_->shell(shell_iter).nfunction()*nocc, naux_, 1.0, met_m_0_5, naux_, u, naux_ * primary_->shell(shell_iter).nfunction(), 0.0, x + naux_ * nbf_ * (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index()), naux_ * primary_->shell(shell_iter).nfunction() );
                }
            }
            break;
        case 'B':
            std::unique_ptr<double[]> PHI(new double[naux_]);
            double* phi = PHI.get();
            double* met_m_1_0 = get_metric_power(-1.0);
// This is the coulomb matrix. We will assume that this is zeroed in the
//     wrapping function
            printf("ping XBLOCK\n");
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
                printf("ping XBLOCK\n");
                compute_dense_AO_block_p_pQq(shell_iter, ao_block, eri);
                printf("ping XBLOCK\n");
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
                    C_DGEMV('N', (int) naux_, (int) nbf_, 1.0, ao_block + func_it*naux_*nbf_, nbf_, d + primary_->shell(shell_iter).function_index()*nbf_ + func_it * nbf_, 1, 1.0, coulomb_vector, 1 );
                    printf("ping XBLOCK\n");
                }
                if (compute_k) {
                    C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, nbf_, 1.0, ao_block, nbf_, c, nocc, 0.0, u, nocc);
                    printf("ping XBLOCK\n");
                    C_DGEMM('N', 'N', naux_, primary_->shell(shell_iter).nfunction()*nocc, naux_, 1.0, met_m_0_5, naux_, u, naux_ * primary_->shell(shell_iter).nfunction(), 0.0, x + naux_ * nbf_ * (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index()), naux_ * primary_->shell(shell_iter).nfunction() );
                    printf("ping XBLOCK\n");
                }
            }
            printf("ping XBLOCK\n");
            C_DGEMV('N', (int) naux_, (int) naux_, 1.0, met_m_1_0, naux_, coulomb_vector, 1, 1.0, phi, 1);
            for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
                compute_dense_AO_block_p_pQq(shell_iter, ao_block, eri);
                for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
                    C_DGEMV('T', (int) naux_, (int) nbf_, 1.0, ao_block + func_it * naux_ * nbf_ , nbf_, phi, 1, 1.0,  j + naux_*nbf_ * (primary_->shell(shell_iter).function_index() + func_it ), 1 );
                }
            }
            C_DCOPY(naux_, phi, 1, coulomb_vector, 1);
            break;
    }
}

void DirectDFJK::unprune_j_pQq( size_t mu, double* pruned_j, double* j){
	for (size_t lambda = 0; lambda < nbf_; lambda++) {
		if (schwarz_func_map_pQq_[mu*nbf_ + lambda]){
			j[mu * nbf_ + lambda] = pruned_j[mu*nbf_ + schwarz_func_map_pQq_[lambda]];
		}
	}
}

} //namespace psi
