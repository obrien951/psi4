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

#include "psi4/libmints/siminteri.h"
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

//buffer[0] is the simint buffer
//buffer[1] is the regular buffer


	printf("auxiliary shells %zu \n", auxiliary_->nshell());
	printf("primary shells %zu \n", primary_->nshell());

	for (int i = 0; i < primary_->nshell(); i++) {
		printf("primary->shell(%d).nfunction() is %zu \n", i, primary_->shell(i).nfunction());
	}

	for (int i = 0; i < auxiliary_->nshell(); i++) {
		printf("auxiliary->shell(%d).nfunction() is %zu \n", i, auxiliary_->shell(i).nfunction());
	}

	std::shared_ptr<BasisSet> zero_ = BasisSet::zero_ao_basis_set();

	auto rifactory = std::make_shared<IntegralFactory>(auxiliary_, zero_, primary_, primary_);
	std::vector<std::shared_ptr<SimintTwoElectronInt>> siminteri(1);

	siminteri[0] = std::make_shared<SimintTwoElectronInt>(rifactory.get(), 0, true);

	printf("aux shell0 functions %zu \n", auxiliary_->shell(0).nfunction());
	printf("pri shell0 functions %zu \n", primary_->shell(0).nfunction());


	std::vector<std::shared_ptr<TwoBodyAOInt>> eri(1);


	std::vector<const double*> buffer(2);
	buffer[0] = siminteri[0]->buffer();

	//siminteri[0]->compute_shell_blocks(0,0,1,1);

	printf("simint calculated AO is %f \n", buffer[0][0]);


	eri[0] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
	buffer[1] = eri[0]->buffer();

	eri[0]->compute_shell(0,0,0,0);

	printf("libmints calculated AO is %f \n", buffer[1][0]);
	printf("Q_shells is %zu\n", Q_shells_);
	printf("p_shells is %zu\n", p_shells_);
	printf("p_shells*p_shells is %zu \n", p_shells_*p_shells_);

	for (int j = 0; j < 1; j++ ) {
		printf("j is %d\n", j);
		siminteri[0]->compute_shell_blocks( 0 , 1  ,1, 28);

		printf("looking at simint integral values\n");
		for (int i = 0; i < 15; i++){
			printf(" %.14f \n", buffer[0][i]);
		}
		printf("\n");
	}

	FILE * eri_ints = fopen("/theoryfs2/ds/obrien/Debug/Psi4/eri_ints.txt","w");

	for (int i = 0; i < primary_->nshell(); i++){
		for (int j = 0; j < primary_->nshell(); j++){
		eri[0]->compute_shell(0,0,i,j);
		fprintf(eri_ints,"i is %d j is %d", i,j);
			for (int inti = 0; inti < primary_->shell(i).nfunction(); inti++) {
				for (int intj = 0; intj < primary_->shell(j).nfunction(); intj++){
					fprintf(eri_ints," %.14f", buffer[1][inti* primary_->shell(j).nfunction() + intj]);
				}
				fprintf(eri_ints,"\n");
			}
		}
	}
	fprintf(eri_ints,"\n");
	fclose(eri_ints);

	for ( int shell_P_it = 0; shell_P_it < auxiliary_->nshell(); shell_P_it++ ) {
		for ( int shell_M_it = 0; shell_M_it < primary_->nshell(); shell_M_it++ ) {
			for ( int shell_N_it = 0; shell_N_it < primary_->nshell(); shell_N_it++ ) {
				
			}
		}
	}



	FILE * tic = fopen("/theoryfs2/ds/obrien/Debug/Psi4/eri_ints.txt","w");



	fclose(tic);

//	int* five = (int*) malloc(sizeof(int));
//	int* twenty_four = (int*) malloc(sizeof(int));
//                                          
//	five[0] = 5;
//	twenty_four[0] = 24;
//
//	double* jjp;
//	double* jkp;
//	double* mjp;
//	double* mkp;
//	double* jdp;
//	double* kdp;
//
//	SharedMatrix Joe__J = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix Joe__K = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix Matt_J = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix Matt_K = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix jdiff = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//	SharedMatrix kdiff = std::make_shared<Matrix>(1, twenty_four, twenty_four, 0);
//
//
//
//
//
//
//
//
//	Matt_J->load("/theoryfs2/ds/obrien/Debug/Psi4/memdfjk_j.txt");
//	Matt_K->load("/theoryfs2/ds/obrien/Debug/Psi4/memdfjk_k.txt");
//
////	SharedMatrix tmp_c = std::make_shared<Matrix>(1, twenty_four, five, 0);
//	
////	tmp_c->load("/theoryfs2/ds/obrien/Debug/Psi4/C_memdfjk.txt");
//
////	C_left_ao_.push_back(tmp_c);
////	C_right_ao_.push_back(tmp_c);
////	compute_D_ao(0);
//
////	std::stringstream J_name;
////	std::stringstream K_name;
////   	J_name << "J " << 0 << " (AO)";
////   	K_name << "K " << 0 << " (AO)";
////	J_ao_.push_back(std::make_shared<Matrix>(J_name.str(), (int) nbf_, (int) nbf_));
////	K_ao_.push_back(std::make_shared<Matrix>(K_name.str(), (int) nbf_, (int) nbf_));
//	
////	num_blocks_ = 2;
//
////	Block_funcs_.resize(0);
//
////	Block_funcs_.push_back(0);
////	Block_funcs_.push_back(0);
//
////	Shell_starts_.push_back(0);
////	Shell_starts_.push_back((Q_shells_/2) + 1);
//
////	Shell_stops_.push_back(Q_shells_/2);
////	Shell_stops_.push_back(Q_shells_ - 1);
//
////	for (size_t i = 0; i <= Q_shells_/2; i++) {
////		Block_funcs_[0] += auxiliary_->shell(i).nfunction();
////	}
//
////	for (size_t i = (Q_shells_/2) + 1; i < Q_shells_; i++) {
////		Block_funcs_[1] += auxiliary_->shell(i).nfunction();
////	}
//
////	biggest_block_ = std::max( Block_funcs_[1], Block_funcs_[0]);
//
////	biggest_block_ = biggest_block_ * nbf_*nbf_;
//
////	compute_JK();
//
//	Joe__J->load("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_J.txt");
//	Joe__K->load("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_K.txt");
//
//	jjp = Joe__J->pointer()[0];
//	jkp = Joe__K->pointer()[0];
//	mjp = Matt_J->pointer()[0];
//	mkp = Matt_K->pointer()[0];
//	jdp = jdiff->pointer()[0];
//	kdp = kdiff->pointer()[0];
//
//
//	for (int i = 0; i < 24; i++){
//		for (int j = 0; j < 24; j++){
//			jdp[i*24+j] = kdp[i*24+j] = 0.0;
//			jdp[i*24+j] = jjp[i*24+j] - mjp[i*24+j];
//			kdp[i*24+j] = jkp[i*24+j] - mkp[i*24+j];
//		}
//	}
//
//	jdiff->save("/theoryfs2/ds/obrien/Debug/Psi4/diff_J.txt", false, false, true);
//	kdiff->save("/theoryfs2/ds/obrien/Debug/Psi4/diff_K.txt", false, false, true);
//
//	free(five);
//	free(twenty_four);
//	sparsity_prep_pQq();

/*
	double memories[300] = {19660800.0, 39321600.0, 58982400.0, 78643200.0, 98304000.0, 117964800.0, 137625600.0, 157286400.0, 176947200.0, 196608000.0, 216268800.0, 235929600.0, 255590400.0, 275251200.0, 294912000.0, 314572800.0, 334233600.0, 353894400.0, 373555200.0, 393216000.0, 412876800.0, 432537600.0, 452198400.0, 471859200.0, 491520000.0, 511180800.0, 530841600.0, 550502400.0, 570163200.0, 589824000.0, 609484800.0, 629145600.0, 648806400.0, 668467200.0, 688128000.0, 707788800.0, 727449600.0, 747110400.0, 766771200.0, 786432000.0, 806092800.0, 825753600.0, 845414400.0, 865075200.0, 884736000.0, 904396800.0, 924057600.0, 943718400.0, 963379200.0, 983040000.0, 1002700800.0, 1022361600.0, 1042022400.0, 1061683200.0, 1081344000.0, 1101004800.0, 1120665600.0, 1140326400.0, 1159987200.0, 1179648000.0, 1199308800.0, 1218969600.0, 1238630400.0, 1258291200.0, 1277952000.0, 1297612800.0, 1317273600.0, 1336934400.0, 1356595200.0, 1376256000.0, 1395916800.0, 1415577600.0, 1435238400.0, 1454899200.0, 1474560000.0, 1494220800.0, 1513881600.0, 1533542400.0, 1553203200.0, 1572864000.0, 1592524800.0, 1612185600.0, 1631846400.0, 1651507200.0, 1671168000.0, 1690828800.0, 1710489600.0, 1730150400.0, 1749811200.0, 1769472000.0, 1789132800.0, 1808793600.0, 1828454400.0, 1848115200.0, 1867776000.0, 1887436800.0, 1907097600.0, 1926758400.0, 1946419200.0, 1966080000.0, 1985740800.0, 2005401600.0, 2025062400.0, 2044723200.0, 2064384000.0, 2084044800.0, 2103705600.0, 2123366400.0, 2143027200.0, 2162688000.0, 2182348800.0, 2202009600.0, 2221670400.0, 2241331200.0, 2260992000.0, 2280652800.0, 2300313600.0, 2319974400.0, 2339635200.0, 2359296000.0, 
	2378956800.0, 2398617600.0, 2418278400.0, 2437939200.0, 2457600000.0, 2477260800.0, 2496921600.0, 2516582400.0, 2536243200.0, 2555904000.0, 2575564800.0, 2595225600.0, 2614886400.0, 2634547200.0, 2654208000.0, 2673868800.0, 2693529600.0, 2713190400.0, 2732851200.0, 2752512000.0, 2772172800.0, 2791833600.0, 2811494400.0, 2831155200.0, 2850816000.0, 2870476800.0, 2890137600.0, 2909798400.0, 2929459200.0, 2949120000.0, 2968780800.0, 2988441600.0, 3008102400.0, 3027763200.0, 3047424000.0, 3067084800.0, 3086745600.0, 3106406400.0, 3126067200.0, 3145728000.0, 3165388800.0, 3185049600.0, 3204710400.0, 3224371200.0, 3244032000.0, 3263692800.0, 3283353600.0, 3303014400.0, 3322675200.0, 3342336000.0, 3361996800.0, 3381657600.0, 3401318400.0, 3420979200.0, 3440640000.0, 3460300800.0, 3479961600.0, 3499622400.0, 3519283200.0, 3538944000.0, 3558604800.0, 3578265600.0, 3597926400.0, 3617587200.0, 3637248000.0, 3656908800.0, 3676569600.0, 3696230400.0, 3715891200.0, 3735552000.0, 3755212800.0, 3774873600.0, 3794534400.0, 3814195200.0, 3833856000.0, 3853516800.0, 3873177600.0, 3892838400.0, 3912499200.0, 3932160000.0, 3951820800.0, 3971481600.0, 3991142400.0, 
	4010803200.0, 4030464000.0, 4050124800.0, 4069785600.0, 4089446400.0, 4109107200.0, 4128768000.0, 4148428800.0, 4168089600.0, 4187750400.0, 4207411200.0, 4227072000.0, 4246732800.0, 4266393600.0, 4286054400.0, 4305715200.0, 4325376000.0, 4345036800.0, 4364697600.0, 4384358400.0, 4404019200.0, 4423680000.0, 4443340800.0, 4463001600.0, 4482662400.0, 4502323200.0, 4521984000.0, 4541644800.0, 4561305600.0, 4580966400.0, 4600627200.0, 4620288000.0, 4639948800.0, 4659609600.0, 4679270400.0, 4698931200.0, 4718592000.0, 4738252800.0, 4757913600.0, 4777574400.0, 4797235200.0, 4816896000.0, 4836556800.0, 4856217600.0, 4875878400.0, 4895539200.0, 4915200000.0, 4934860800.0, 4954521600.0, 4974182400.0, 4993843200.0, 5013504000.0, 5033164800.0, 5052825600.0, 5072486400.0, 5092147200.0, 5111808000.0, 5131468800.0, 5151129600.0, 5170790400.0, 5190451200.0, 5210112000.0, 5229772800.0, 5249433600.0, 5269094400.0, 5288755200.0, 5308416000.0, 5328076800.0, 5347737600.0, 5367398400.0, 5387059200.0, 5406720000.0, 5426380800.0, 5446041600.0, 5465702400.0, 5485363200.0, 5505024000.0, 5524684800.0, 5544345600.0, 5564006400.0, 5583667200.0, 5603328000.0, 5622988800.0, 5642649600.0, 5662310400.0, 5681971200.0, 5701632000.0, 5721292800.0, 5740953600.0, 5760614400.0, 5780275200.0, 5799936000.0, 5819596800.0, 5839257600.0, 5858918400.0, 5878579200.0, 5898240000.0};

	for (int i = 0; i < 300; i++){
		memory_ = static_cast<size_t>(memories[i]);
		prepare_p_blocks();
        printf("memory_ is %zu Mib num_blocks_ is %zu\n", memory_ * 4 * 8L / (3* 1024L * 1024L) ,num_blocks_ );
	}
*/	
		

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
for (int i = 0; i < 1000; i ++) c[i] = 0.0;

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

	C_DGEMM('N', 'N', 4, 5, 6, 1.0, A, 6, B, 5, 0.0, C, 10);

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
        outfile->Printf("    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
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
	sparsity_prep_pQq();

	//prepare blocks tells us at which indices we start and stop each block.
	// I'm making it a separate function because I've seen evidence that 
	// DFHelper does some resizing, and I want to be able to match that.
//dirty
	biggest_shell_ = 0;
	for (int i = 0; i < primary_->nshell(); i++){
		if (primary_->shell(i).nfunction() > biggest_shell_) {
			biggest_shell_ = primary_->shell(i).nfunction();
		}
	}
//end dirty

	// prepares the coulomb metric. Of course, our calculations depend on the
	//  inverse of the coulomb metric and the (-0.5) power of the metric, so
	//  we will prepare those immediately afterwards
	//prepare_metric();
	get_met();
	
	prepare_metric_power(-0.5);

	//for what it's worth, we should have the density matrix
	//prepare_D_symm();

}

void sparsity_prep_Qpq() { }

void DirectDFJK::sparsity_prep_pQq(){
	printf("started sparsity prep\n");
	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
#endif

	// variables to hold data
	double global_max_int = 0.0;


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
    int rank = 0;
#ifdef _OPENMP
    rank = omp_get_thread_num();
#endif
    eri[rank] = std::shared_ptr<TwoBodyAOInt>(rifactory->eri());
	buffer[rank] = eri[rank]->buffer();
    }

	// Sparsity member data
	schwarz_func_ints_.resize(nbf_);

	std::vector<size_t> dense_its(nbf_);

	//AO shell metadata
	size_t Shell_M_count;
	size_t Shell_M_start;
	size_t Shell_N_count;
	size_t Shell_N_start;
	double val = 0.0;

#pragma omp parallel for schedule(guided) num_threads(procs) private( Shell_M_count, Shell_M_start, Shell_N_count, Shell_N_start, val) reduction(max:global_max_int)
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
//					integrals[ (Shell_N_start + func_n) * nbf_ + Shell_M_start + func_m ]
//					=
					val 
					= 
 fabs(buffer[rank][ func_m * Shell_N_count * Shell_M_count * Shell_N_count + 
				    func_n * Shell_M_count * Shell_N_count + 
				    func_m * Shell_N_count + 
				    func_n]);
					if (global_max_int < val) global_max_int = val;
					
					if (shel_maxes[Shell_M * p_shells_ + Shell_N] < val) {
						shel_maxes[Shell_M * p_shells_ + Shell_N]
//						=
//						shel_maxes[Shell_N * p_shells_ + Shell_M]
						= val;
					}
				}
			}
		}
	}

	tolerance_ = cutoff_ * cutoff_ / global_max_int;
	schwarz_shell_mask_pQq_.resize(p_shells_);
	schwarz_func_starts_pQq_.resize(p_shells_);
	schwarz_dense_funcs_.resize(p_shells_);
	size_t count;

	for ( size_t shell_it_out = 0U; shell_it_out < p_shells_; shell_it_out++ ) {
		count = 0U;
		for (size_t shell_it_in = 0U; shell_it_in < p_shells_; shell_it_in++ ) {
			if(shel_maxes[shell_it_out*p_shells_ + shell_it_in] >  tolerance_) {
				schwarz_shell_mask_pQq_[shell_it_out].push_back(shell_it_in);
				schwarz_func_starts_pQq_[shell_it_out].push_back(count);
				count += primary_->shell(shell_it_in).nfunction();
			}
		}
		schwarz_dense_funcs_[shell_it_out] = count;
	}
    printf("finished sparsity prep\n");
}

// calculates a block of AO's then conducts the proper contractions to get them
//   where they need to go. I don't think it makes sense to have seperate J and
//   K function calls as J and K are never fully built in a single function.
void DirectDFJK::compute_JK() {
	printf("starting compute_JK \n");
	BB_ = false;
//printf("num blocks is %zu\n", num_blocks_);
//I'm leaving this commented out so that nobody tries to write this code back
//  in. the number of blocks is caluclated with each hartree fock iteration.
//  it's a short code that does it, so it doesn't need to be eliminated. However
//  , it is necessary to have the correct block size for the SAD guess.
	if ( pQq_ ) {
//    	pQp();
			//pQp_sparse();

		if (nocc_last_ == C_left_ao_[0]->ncol() ) {
			printf("about to start pQp_mn_mP_sparse()\n");
			outfile->Printf("we're going to use mP sparsity as set\n");
			pQp_mn_mP_sparse();
			printf("just finished pQp_mn_mP_sparse()\n");
		} else {
			printf("about to start pQp_mn_sparse_set_mP()\n");
			outfile->Printf("we're going to set mP sparsity\n");
			pQp_mn_sparse_set_mP();
			printf("just finished pQp_mn_sparse_set_mP()\n");
		}

	}  else {

		if (num_blocks_ == 1) {
			build_jk_CC_Qpq_direct();
		} else {
    		build_jk_CC_Qpq_blocks();
		}

	}
	printf("finishing compute_JK\n");
}

void DirectDFJK::postiterations() {
}

// Determines the maximum amount of memory that DirectDFJK will be using
//   at any one time. To be sure, the memory used to construct AO blocks
//   will only be allocated during the compute_JK() method, so it will
//   be available during other procedures. 
// Uses blocking and sparsity information to set total_needs which is 
//   returned by memory_estimate.
void DirectDFJK::our_needs(){
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
	if (pQq_) { 
//		sparsity_prep_pQq();
//		prepare_p_blocks();
	}

	if (Qpq_) {
		if (Qpq_store_sparse_){ sparsity_prep_Qpq(); }
		prepare_Q_blocks();
	}

	our_needs();
	return total_needs_;
}

void DirectDFJK::prepare_metric_power(double power){
	SharedMatrix my_met;
	auto coul_met = std::make_shared<FittingMetric>(auxiliary_, true);
	coul_met->form_fitting_metric();
	my_met = coul_met->get_metric();
	if ( !(1e-13 > fabs(1.0 - power)) ){
		my_met->power(power, condition_);
	}
	
	met_cols_.push_back(my_met->ncol());
	met_rows_.push_back(my_met->nrow());
	metric_.push_back(my_met);
	met_powers_.push_back(power);

    while (power < 0.0) {
        power = power + 1.0;
    }
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

// Computes the LU decomposition of the coulomb metric and stores it in
//   member vector CMPQ_LU_.
//   pivots are storec in vector PERMUTE_
void DirectDFJK::get_met(){
	auto coul_met = std::make_shared<FittingMetric>(auxiliary_, true);

	coul_met->form_fitting_metric();
	SharedMatrix my_met;
	my_met = coul_met->get_metric();
	PERMUTE_.resize(naux_);
	CMPQ_LU_.resize(naux_*naux_);

	double* cmpq_lup = &CMPQ_LU_.front();
	double* metp = my_met->pointer()[0];
	int* pert = &PERMUTE_.front();
	
	C_DGETRF( naux_, naux_, metp, naux_, pert);
	C_DCOPY(naux_*naux_, metp, 1, cmpq_lup, 1);
	
}

// called after sparsity_prep_pQq() which is called in the memory estimate
// is called in the memory estimate because the maximum permissable 
// block size is needed for an effective memory estimate, and determining
// that is equivalent to estimating the total amount of memory
// Determines the number of blocks and 
void DirectDFJK::prepare_p_blocks() {
	Block_funcs_.resize(0);
	Shell_starts_.resize(0);
	Shell_stops_.resize(0);
	k_disps_.resize(0);

	size_t current_costs = 0;
	// block max functions
	size_t bmf = 0;

	size_t block_size = 0;
	size_t shell_funcs;
    size_t biggest_shell = 0;

	// one of two things. Either nelectron/2 + 1 or C_left_ao->ncol()
	size_t charges_z;

	if (C_left_ao_.size() == 0) {
	// atomic charges float
	double charges_f;

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
	} else {
		charges_z = C_left_ao_[0]->ncol();
	}
//
// D_ao_ and D_
	current_costs += 3*nbf_*nbf_;
// Coulomb Metrics
	current_costs += 2*naux_*naux_;
	current_costs += 4*nbf_*nbf_;
// adding in costs for U and A
    current_costs += biggest_shell*naux_*charges_z;
    current_costs += biggest_shell*naux_*charges_z;

//	if ( memory_ - current_costs > memory_) {
//		std::stringstream error_message;
//		error_message << "DirectDFJK: not enough memory for overhead! Overhead requires " << current_costs * 8 / (1024*1024) << " [MiB].";
//		throw PSIEXCEPTION(error_message.str().c_str());
//	}
   
	Shell_starts_.push_back(0);
    
	for (size_t shell_iter = 0; shell_iter < p_shells_; shell_iter++) {
		shell_funcs = primary_->shell(shell_iter).nfunction();
// The 5 is where we're getting our fudge factor we're getting charges_z
//   to stand in for n_occ, and there should only be a factor of 2 included for
//   the number of tensors with dimension occ. The extra 3 is a massive amount of 
//   memory that should cover the overhead from the reast of the calculation.
//   charges_z also includes a fudge factor
		if (  (block_size + shell_funcs + 2) * naux_ * ( 2*charges_z )  >  memory_ - current_costs ){
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
	Shell_stops_.push_back(p_shells_ - 1);

	Block_funcs_.push_back(block_size);

	biggest_block_ = bmf * naux_ * nbf_;
	num_blocks_ = Block_funcs_.size();
    k_disps_.resize(num_blocks_);
	size_t col_disp = 0U;
	size_t row_disp = 0U;

	for (size_t row_it = 0; row_it < num_blocks_; row_it++) {
		k_disps_[row_it].resize(num_blocks_);
		for (size_t col_it = 0; col_it < num_blocks_; col_it++) {
			k_disps_[row_it][col_it] = nbf_* row_disp + col_disp;
			col_disp += Block_funcs_[col_it];
		}
		row_disp += Block_funcs_[row_it];
		col_disp = 0U;
	}

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
#ifdef _OPENMP
	procs = omp_nthread_;
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
    
	int rank = 0;
    size_t shell_count = primary_->shell(shell).nfunction();
    size_t ShellP_count;
	size_t ShellP_start;
    size_t ShellN_count;
	size_t ShellN_start; 
    size_t Shell_ind_0;
    size_t Buff_ind_0;
// the value put in mat_size is nbf*naux
    size_t mat_size = nbf_ * naux_;
    
// In Matt's Code and Rob's code (and, indeed, some of mine), we see 3 loops.
//   however, in this code, we only need two because we only have one shell
//   to compute in the slowest running index, so we're just going to write
//   2 loops to accomodate that.
#pragma omp parallel for schedule(guided) num_threads(procs) private(ShellP_count, ShellP_start, ShellN_count, ShellN_start, Shell_ind_0, Buff_ind_0, rank)
    for (size_t ShellP = 0; ShellP < Q_shells_; ShellP++) {

#ifdef _OPENMP
		rank = omp_get_thread_num();
#endif
        ShellP_count = auxiliary_->shell(ShellP).nfunction();
		ShellP_start = auxiliary_->shell(ShellP).function_index();
        for (size_t ShellN = 0; ShellN < p_shells_; ShellN++ ) {
            ShellN_count = primary_->shell(ShellN).nfunction();
			ShellN_start = primary_->shell(ShellN).function_index();
            Buff_ind_0 = ShellN_count * shell_count;
			eri[rank]->compute_shell(ShellP, 0, shell, ShellN);
            int i = 0;
            for (size_t func_m  = 0; func_m < shell_count; func_m++){
                for (size_t func_p = 0; func_p < ShellP_count; func_p++ ) {
                    for (size_t func_n = 0; func_n < ShellN_count; func_n++) {
                        ao_block[func_m * nbf_*naux_+ (ShellP_start + func_p) * nbf_+ ShellN_start + func_n]
                         =
                        buffer[rank][func_p * Buff_ind_0 + func_m * ShellN_count + func_n];
                    }
                }
            }
        }
    }
}

void DirectDFJK::compute_sparse_AO_block_p_pQq( size_t shell, double* ao_block, std::vector<std::shared_ptr<TwoBodyAOInt>> eri){
	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
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
	size_t Shell_count = primary_->shell(shell).nfunction();
	size_t ShellN_start;
	size_t ShellN_count;
	size_t Shell_ind_0;
	size_t func_starts_sum;// schwarz_func_starts_pQq_

// If we try to save something at &ao_block + ShellP_start *nbv_^2, we
//   will get a segfault in all but the 0th block. disp short for displacemtnt


// We're building A_\mu P \nu

// Loop over auxes
#pragma omp parallel for schedule(guided) num_threads(procs) private(ShellP_start, ShellP_count, ShellN_start, ShellN_count, Shell_ind_0)
	for ( size_t ShellP = 0; ShellP < Q_shells_; ShellP++ ) {
		int rank = 0;
#ifdef _OPENMP
		rank = omp_get_thread_num();
#endif
    	ShellP_start = auxiliary_->shell(ShellP).function_index();
    	ShellP_count = auxiliary_->shell(ShellP).nfunction();
// Inner Loop over prims
		for ( size_t ShellN = 0; ShellN < schwarz_shell_mask_pQq_[shell].size() ; ShellN++ ) {
			ShellN_start = schwarz_func_starts_pQq_[shell][ShellN]; 
// primary_->shell(schwarz_shell_mask_pQq_[shell][ShellN]).function_index();
			ShellN_count = primary_->shell(schwarz_shell_mask_pQq_[shell][ShellN]).nfunction();
			Shell_ind_0 = Shell_count*ShellN_count;
// Calculate AO's
			eri[rank]->compute_shell(ShellP, 0, shell, schwarz_shell_mask_pQq_[shell][ShellN]);
// repack AO's
			for( size_t intp = 0; intp < ShellP_count; intp++ ){
				for( size_t intm = 0; intm < Shell_count; intm++ ){
					for( size_t intn = 0; intn < ShellN_count; intn++ ){
	ao_block[ intm * schwarz_dense_funcs_[shell] * naux_ + (intp + ShellP_start) * schwarz_dense_funcs_[shell] +  ShellN_start+ intn ] 
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
void DirectDFJK::compute_AO_block_Qpq( size_t start_Q, size_t stop_Q, double* ao_block, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ){
	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
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

void DirectDFJK::compute_AO_block_p_pQq_mn_sparse_set_mP( size_t shell, double* ao_block, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ) {
	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
#endif

	std::unique_ptr<double[]> T_AO(new double[biggest_shell_* nbf_ * naux_]);
	double * temp_ao = T_AO.get();

	std::vector<const double*> buffer(procs);
    std::vector<std::vector<bool>> colc(procs); 
	std::vector<bool> shell_track(procs);

#pragma omp parallel num_threads(procs)
	{
	int rank = 0;
#ifdef _OPENMP
	rank = omp_get_thread_num();
#endif
	buffer[rank] = eri[rank]->buffer();
    colc[rank].resize(biggest_shell_);
	}

 // ERI indexing Variables
 // Incides are done with romanization of greek letters \Mu -> M \Nu ->N
  	size_t ShellP_start;
    size_t ShellP_count;
 	size_t ShellM_start = primary_->shell(shell).function_index();
 	size_t ShellM_count = primary_->shell(shell).nfunction();
 	size_t ShellN_start;
 	size_t ShellN_count;
 	size_t Shell_ind_0;

	printf("done initiating variables\n");
	printf("about to calculate aos\n");

	printf("%zu\n", Q_shells_);

	for ( size_t ShellP = 0; ShellP < Q_shells_; ShellP++ ) {
		printf("ShellP is %zu\n", ShellP);

    	int rank = 0;
#ifdef _OPENMP
    	rank = omp_get_thread_num();
#endif
		ShellP_start = auxiliary_->shell(ShellP).function_index();
		ShellP_count = auxiliary_->shell(ShellP).nfunction();
//initialize terms used for mP sparsity
//using a for loop as opposed to a parallel block because the outer 
//        loop is parallelized
        for (size_t i = 0; i < ShellP_count; i++) { colc[rank][i] = false; }

		for ( size_t ShellN = 0; ShellN < schwarz_shell_mask_pQq_[shell].size(); ShellN++ ) {
			ShellN_start = primary_->shell(ShellN).function_index();
			ShellN_count = primary_->shell(ShellN).nfunction();
            Shell_ind_0 = ShellN_count * ShellM_count;
			eri[rank]->compute_shell(ShellP, 0, shell, schwarz_shell_mask_pQq_[shell][ShellN]);
				for ( size_t intm = 0; intm < ShellM_count; intm++ ) {
			for ( size_t intp = 0; intp < ShellP_count; intp++ ) {
					for ( size_t intn = 0; intn < ShellN_count; intn++ ) {
						temp_ao[ intm*naux_*schwarz_shell_mask_pQq_[shell].size() 
                                + (intp + ShellP_start)*schwarz_shell_mask_pQq_[shell].size()
                                + intn  ]
                        =
                        buffer[rank][ intp * Shell_ind_0
                                    + intm * ShellN_count
                                    + intn ];

                        if ( buffer[rank][ intp * Shell_ind_0 + intm * ShellN_count + intn ] > tolerance_ ) { colc[rank][intp] = true; }
					}
				}
			}
		}
// will need to sort the sparsity maps for functiosn and shells
//     for tensor contractions
        for ( size_t i = 0; i < ShellP_count; i++ ) { 
            if (colc[rank][i]) { 
                mP_func_map_pQq_[shell].push_back(ShellP_start + i);
               	shell_track[rank] = true;
            }
        }
		if ( shell_track[rank] ) {
			mP_shel_map_pQq_[shell].push_back(ShellP);
		}
	}

	printf("done calculating aos\n");

    std::sort(mP_func_map_pQq_[shell].begin(), mP_func_map_pQq_[shell].end(), DirectDFJK::sztcmp );
    std::sort(mP_shel_map_pQq_[shell].begin(), mP_shel_map_pQq_[shell].end(), DirectDFJK::sztcmp );

	printf("done sorting\n");

	size_t mP_size = mP_func_map_pQq_[shell].size();
	size_t mn_size = schwarz_dense_funcs_[shell];
	size_t nu_base;

	printf("done setting pruning variables\n");

for (int i = 0; i < mP_size; i++) {
		printf("%zu \n",mP_func_map_pQq_[shell][i] );
	}



// We really want to get our clock-cycles' worth from
//    this sparsity. we're going to repack the AOs for a sparse
//    tensor contraction.
#pragma omp parallel for num_threads(procs)
	for (size_t mu_it = 0; mu_it < primary_->shell(shell).nfunction(); mu_it++ ) {
		for ( size_t Q_it = 0; Q_it < mP_size; Q_it++ ) {
			for (size_t nu_it = 0; nu_it < mn_size; nu_it++ ) {
					ao_block[ mu_it * mP_size * mn_size
							+ Q_it  * mn_size
							+ nu_it ]
							=
					temp_ao [ mu_it * naux_ * mn_size
							+ Q_it * mn_size 
							+ nu_it];
			}
		}
	}
}

void DirectDFJK::compute_AO_block_p_pQq_mn_mP_sparse( size_t shell, double* ao_block, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ) {
	printf("in compute_AO_block_p_pQq_mn_mP_sparse\n");

	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
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
 	size_t ShellM_start = primary_->shell(shell).function_index();
 	size_t ShellM_count = primary_->shell(shell).nfunction();
 	size_t ShellN_start;
 	size_t ShellN_count;
 	size_t Shell_ind_0; 


	for ( size_t ShellP = 0; ShellP < mP_shel_map_pQq_[shell].size(); ShellP++ ) {
    	int rank = 0;
#ifdef _OPENMP
    	rank = omp_get_thread_num();
#endif
		ShellP_start = auxiliary_->shell(mP_shel_map_pQq_[shell][ShellP]).function_index();
		ShellP_count = auxiliary_->shell(mP_shel_map_pQq_[shell][ShellP]).nfunction();
//initialize terms used for mP sparsity
//using a for loop as opposed to a parallel block because the outer 
//        loop is parallelized

		for ( size_t ShellN = 0; ShellN < schwarz_shell_mask_pQq_[shell].size(); ShellN++ ) {
			ShellN_start = primary_->shell(ShellN).function_index();
			ShellN_count = primary_->shell(ShellN).nfunction();
            Shell_ind_0 = ShellN_count * ShellM_count;
eri[rank]->compute_shell( mP_shel_map_pQq_[shell][ShellP], 0, shell, schwarz_shell_mask_pQq_[shell][ShellN]);
				for ( size_t intm = 0; intm < ShellM_count; intm++ ) {
			for ( size_t intp = 0; intp < ShellP_count; intp++ ) {
					for ( size_t intn = 0; intn < ShellN_count; intn++ ) {
						ao_block[ intm*naux_*schwarz_shell_mask_pQq_[shell].size() 
                                + (intp + ShellP_start)*schwarz_shell_mask_pQq_[shell].size()
                                + intn  ]
                        =
                        buffer[rank][ intp * Shell_ind_0
                                    + intm * ShellN_count
                                    + intn ];

					}
				}
			}
		}
	}
	printf("exiting compute_AO_block_p_pQq_mn_mP_sparse\n");
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
	procs = omp_nthread_;
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


	std::unique_ptr<double[]> A(new double[nbf_*nbf_*naux_]);
	std::unique_ptr<double[]> B(new double[nbf_*nbf_*naux_]);
	std::unique_ptr<double[]> U(new double[x_size_]);
	std::unique_ptr<double[]> X(new double[x_size_]);
	std::unique_ptr<double[]> V(new double[naux_]);
	std::unique_ptr<double[]> PHI(new double[naux_]);

	double* a = A.get();
	double* b = B.get();
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

timer_on("Tensor Contractions");
	C_DGEMM( 'N', 'N', naux_, nbf_*nbf_, naux_, 1.0, met_m_0_5, naux_, a, nbf_*nbf_, 0.0, b, nbf_*nbf_);

	C_DGEMV( 'N', naux_, nbf_*nbf_, 1.0, b, nbf_*nbf_, d, 1, 0.0, v, 1 );

	C_DGEMV( 'T', naux_, nbf_*nbf_, 1.0, b, nbf_*nbf_, v, 1, 0.0, j, 1 );

//	V_gets_AD( naux_, v + in_block_off, a, d);
//	C_DGEMV( 'N', (int) naux_, (int) nbf_squared, 1.0,         a, (int) nbf_squared,   d, 1, 0.0,   v, 1);

	// PHI /gets \cmpq V
//	C_DGEMV( 'T', (int) naux_, (int)       naux_, 1.0, met_m_1_0, (int)       naux_,   v, 1, 0.0, phi, 1);

//	Accumulate_J( Block_funcs_[0] , j, a, phi);
// Accumulate_J
//	C_DGEMV( 'T', (int) naux_, (int) nbf_squared, 1.0,         a, (int) nbf_squared, phi, 1, 0.0,   j, 1);

//	U_gets_AC(Block_funcs_[0], static_cast<size_t>(C_left_ao_[0]->ncol()), u, a, c);
// Form U by flattening out U over the auxiliary basis.

//changed a to b. changed a to b. changed a to b. changed a to b.
	C_DGEMM( 'N', 'N', (int) (nbf_*naux_), C_left_ao_[0]->ncol(), (int) nbf_, 1.0, b, (int) nbf_, c, C_left_ao_[0]->ncol(), 0.0, x, C_left_ao_[0]->ncol());
//changed a to b. changed a to b. changed a to b. changed a to b.

//	X_accumulates_JU(Block_funcs_[0], Block_funcs_[0], x_slice_[0], ou_block_off, in_block_off, static_cast<size_t>(C_left_ao_[0]->ncol()), x, met_m_0_5, u );
// Form x by flattening out U along the auxiliary basis and transposing
// to contract over the auxiliary basis set.
//	C_DGEMM( 'T', 'T', ((int) nbf_) * C_left_ao_[0]->ncol(), (int) naux_, (int) naux_, 1.0, u, ((int) nbf_) * C_left_ao_[0]->ncol(), met_m_0_5, (int) naux_, 0.0, x, (int) naux_);

//	Accumulate_K_c_is_c(Block_funcs_[0], x_slice_[0], static_cast<size_t>(C_left_ao_[0]->ncol()), k, x);
//	C_DGEMM( 'N', 'T', (int) nbf_, (int) nbf_, ((int) naux_) * C_left_ao_[0]->ncol(), 1.0, x, ((int) naux_) * C_left_ao_[0]->ncol(), x, ((int) naux_) * C_left_ao_[0]->ncol()  , 0.0, k, (int) nbf_);

	for (size_t R_it = 0; R_it < naux_; R_it++) {
		C_DGEMM('N', 'T', nbf_, nbf_, C_left_ao_[0]->ncol(), 1.0, x+R_it*nbf_*C_left_ao_[0]->ncol(), C_left_ao_[0]->ncol(), x+R_it*nbf_*C_left_ao_[0]->ncol(), C_left_ao_[0]->ncol(), 1.0, k, nbf_);
	}

	//	J_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_J.txt", false, false, true);
	//	K_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_K.txt", false, false, true);

	//	C_left_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/C_directdfjk.txt", false, false, true);
timer_off("Tensor Contractions");
		x_slice_.erase(x_slice_.begin(), x_slice_.end() );
	}


void DirectDFJK::build_jk_CC_Qpq_blocks() {
	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
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
			compute_AO_block_Qpq( Shell_starts_[block_iter_in], Shell_stops_[block_iter_in], a, eri);
			if ( block_iter_ou == 0 ) {
//				V_gets_AD( Block_funcs_[block_iter_in], v + in_block_off, a, d);// + v_phi_add);
				C_DGEMV( 'N', Block_funcs_[block_iter_in], nbf_squared, 1.0, a, nbf_squared,   d, 1, 0.0, v + in_block_off, 1);
			}
			if ( block_iter_ou == 1 ) {
//				Accumulate_J( Block_funcs_[block_iter_in] , j, a, phi + in_block_off);
				C_DGEMV( 'T', Block_funcs_[block_iter_in], nbf_squared, 1.0, a, nbf_squared, phi + in_block_off, 1, 1.0, j, 1);
			}
timer_on("Tensor Contractions");
//			U_gets_AC(Block_funcs_[block_iter_in], static_cast<size_t>(C_left_ao_[0]->ncol()), u, a, c);
			C_DGEMM( 'N', 'N', Block_funcs_[block_iter_in] * nbf_, C_left_ao_[0]->ncol(), nbf_, 1.0, a, nbf_, c, C_left_ao_[0]->ncol(), 0.0, u, C_left_ao_[0]->ncol());
//			X_accumulates_JU(Block_funcs_[block_iter_ou], Block_funcs_[block_iter_in], x_slice_[0], ou_block_off, in_block_off, static_cast<size_t>(C_left_ao_[0]->ncol()), x, met_m_0_5, u );
			C_DGEMM( 'T', 'N', nbf_* C_left_ao_[0]->ncol(), Block_funcs_[block_iter_ou], Block_funcs_[block_iter_in], 1.0, u, nbf_* C_left_ao_[0]->ncol(), met_m_0_5 + (in_block_off * naux_) + ou_block_off, naux_, 1.0, x, Block_funcs_[block_iter_ou]);
			in_block_off += Block_funcs_[block_iter_in];
timer_off("Tensor Contractions");
		}
timer_on("Tensor Contractions");
		C_DGEMM( 'N', 'T', nbf_, nbf_, Block_funcs_[block_iter_ou] * C_left_ao_[0]->ncol(), 1.0, x, Block_funcs_[block_iter_ou] * C_left_ao_[0]->ncol(), x, Block_funcs_[block_iter_ou] * C_left_ao_[0]->ncol(), 1.0, k, nbf_ );
timer_off("Tensor Contractions");
//		Accumulate_K_c_is_c(Block_funcs_[block_iter_ou], x_slice_[0], static_cast<size_t>(C_left_ao_[0]->ncol()), k, x);
		ou_block_off += Block_funcs_[block_iter_ou];
//		x_block = x_block + Block_funcs[block_iter_ou] * x_slice_[0];
	}
//J_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_J.txt", false, false, true);
//K_ao_[0]->save("/theoryfs2/ds/obrien/Debug/Psi4/directdfjk_K.txt", false, false, true);
	x_slice_.erase( x_slice_.begin(), x_slice_.end() ); 
}

void DirectDFJK::pQp(){
// In principle, this function call should be in preiterations or in
//   memory estimator. However, it depends on knowing the memory_
//   value from input which precludes its calling in either of those
//   places.
	prepare_p_blocks();

	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
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

	int nocc = C_left_ao_[0]->ncol();

	double* j = J_ao_[0]->pointer(0)[0];
	double* k = K_ao_[0]->pointer(0)[0];
	
	std::unique_ptr<double[]> A(new double[biggest_shell_*naux_*nbf_]);
	std::unique_ptr<double[]> U(new double[1U]);
	if (BB_) {
		U.reset(new double[biggest_shell_*naux_*nbf_]);
	} else {
		U.reset(new double[biggest_shell_*naux_*nocc]);
	}
	std::unique_ptr<double[]> XN(new double[biggest_block_/nbf_*nocc]);
	std::unique_ptr<double[]> XO(new double[biggest_block_/nbf_*nocc]);
	std::unique_ptr<double[]> V(new double[naux_]);


	double* a = A.get();
	double* u = U.get();
	double* xn = XN.get();
	double* xo = XO.get();
	double* xh;
	double* v = V.get();
	double Zero = 0.0;
	double* zero = &Zero;

	C_DCOPY( nbf_*nbf_, zero, 0, j, 1);
    C_DCOPY( nbf_*nbf_, zero, 0, k, 1);

	C_DCOPY(naux_, zero, 0, v, 1);
	//for (size_t i = 0; i < naux_; i++) { v[i] = 0.0;}
	char first_char = ( num_blocks_ == 1 ? 'B' : 'V');
	

	size_t xo_ind = 0U;
	size_t xn_ind;

	X_Block(first_char, true, 0, a, xo, u, v, eri);
	C_DGEMM('N', 'T', Block_funcs_[0], Block_funcs_[0], naux_*nocc, 1.0, xo, naux_*nocc, xo, naux_*nocc, 1.0, k, nbf_);
	for (size_t xo_iter = 0; xo_iter < num_blocks_ - 1 ; xo_iter++){
		for (size_t xn_iter = 1; xn_iter < num_blocks_ - xo_iter; xn_iter++){
			if ( xo_iter == 0 && xn_iter != num_blocks_ - 1) {
				X_Block('V', true, xn_iter, a, xn, u, v, eri);
				xn_ind = xn_iter;
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 1.0, k + k_disps_[xn_ind][xn_ind], nbf_ );
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k + k_disps_[xo_ind][xn_ind], nbf_ );
			} else if (xo_iter == 0 && xn_iter == num_blocks_ - 1) {
				X_Block('B', true, xn_iter, a, xn, u, v, eri);
				xn_ind = xn_iter;
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 1.0, k + k_disps_[xn_ind][xn_ind], nbf_);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k + k_disps_[xo_ind][xn_ind], nbf_ );
			} else if ( xo_iter == 1 ) {
				X_Block('P', true, xn_iter, a, xn, u, v, eri);
				xn_ind = xn_iter;
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xo_ind], naux_*nocc, 1.0, xn, naux_*nocc, xo, naux_*nocc, 1.0, k + k_disps_[xn_ind][xo_ind], nbf_);
			} else {
				X_Block('N', true, xn_iter, a, xn, u, nullptr, eri);
				xn_ind = xn_iter;
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xo_ind], naux_*nocc, 1.0, xn, naux_*nocc, xo, naux_*nocc, 0.0, k + k_disps_[xn_ind][xo_ind], nbf_);
			}
		}
		xh = xn;
		xn = xo;//XO.get();
		xo = xh;//XN.get();
		xo_ind = xn_ind;
	}
	for (size_t kf_i = 0; kf_i < nbf_; kf_i++){
		for (size_t kf_j = 0; kf_j < kf_i; kf_j++){
			k[kf_i * nbf_ + kf_j] = k[ kf_j *nbf_ + kf_i];
		}
	}
	
	if (first_char != 'B') {
		X_Block('P', false, 0, a, nullptr, u, v, eri);
	}
}

void DirectDFJK::pQp_sparse(){
// In principle, this function call should be in preiterations or in
//   memory estimator. However, it depends on knowing the memory_
//   value from input which precludes its calling in either of those
//   places.
	prepare_p_blocks();

	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
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

	int nocc = C_left_ao_[0]->ncol();

	double* j = J_ao_[0]->pointer(0)[0];
	double* k = K_ao_[0]->pointer(0)[0];
	
	std::unique_ptr<double[]> A(new double[ biggest_shell_ * naux_ * nbf_ ]);
	std::unique_ptr<double[]> U(new double[biggest_shell_ * naux_*nocc]);
	std::unique_ptr<double[]> XN(new double[biggest_block_/nbf_*nocc]);
	std::unique_ptr<double[]> XO(new double[biggest_block_/nbf_*nocc]);
	std::unique_ptr<double[]> V(new double[naux_]);

	double* a = A.get();
	double* u = U.get();
	double* xn = XN.get();
	double* xo = XO.get();
	double* xh;
	double* v = V.get();
	double Zero = 0.0;
	double* zero = &Zero;


	std::unique_ptr<double[]> P_C(new double[nbf_*nocc]);
	std::unique_ptr<double[]> P_D(new double[nbf_]);
	std::unique_ptr<double[]> P_J(new double[nbf_]);


	double* pruned_c = P_C.get();
	double* pruned_d = P_D.get();
	double* pruned_j = P_J.get();


	C_DCOPY( nbf_*nbf_, zero, 0, j, 1);
    C_DCOPY( nbf_*nbf_, zero, 0, k, 1);
	C_DCOPY(naux_, zero, 0, v, 1);
	//for (size_t i = 0; i < naux_; i++) { v[i] = 0.0;}
	char first_char = ( num_blocks_ == 1 ? 'B' : 'V');

	size_t xo_ind = 0U;
	size_t xn_ind;
	
	outfile->Printf("num_blocks_ is %zu\n", num_blocks_);

	X_Block_sparse(first_char, true, 0, pruned_c, pruned_d, a, xo, u, v, pruned_j, eri);
timer_on("DDF pQq big K DGEMM");
	C_DGEMM('N', 'T', Block_funcs_[0], Block_funcs_[0], naux_*nocc, 1.0, xo, naux_*nocc, xo, naux_*nocc, 1.0, k, nbf_);
timer_off("DDF pQq big K DGEMM");
	for (size_t xo_iter = 0; xo_iter < num_blocks_ - 1 ; xo_iter++){
		for (size_t xn_iter = 1; xn_iter < num_blocks_ - xo_iter; xn_iter++){
			if ( xo_iter == 0 && xn_iter != num_blocks_ - 1) {
				X_Block_sparse('V', true, xn_iter, pruned_c, pruned_d, a, xn, u, v, pruned_j, eri);
				xn_ind = xn_iter;
timer_on("DDF pQq big K DGEMM");
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 1.0, k + k_disps_[xn_ind][xn_ind], nbf_ );
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k + k_disps_[xo_ind][xn_ind], nbf_ );
timer_off("DDF pQq big K DGEMM");
			} else if (xo_iter == 0 && xn_iter == num_blocks_ - 1) {
				X_Block_sparse('B', true, xn_iter, pruned_c, pruned_d, a, xn, u, v, pruned_j, eri);
				xn_ind = xn_iter;
timer_on("DDF pQq big K DGEMM");
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 1.0, k + k_disps_[xn_ind][xn_ind], nbf_);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k + k_disps_[xo_ind][xn_ind], nbf_ );
timer_off("DDF pQq big K DGEMM");
			} else if ( xo_iter == 1 ) {
				X_Block_sparse('P', true, xn_iter, pruned_c, pruned_d, a, xn, u, v, pruned_j, eri);
				xn_ind = xn_iter;
timer_on("DDF pQq big K DGEMM");
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xo_ind], naux_*nocc, 1.0, xn, naux_*nocc, xo, naux_*nocc, 1.0, k + k_disps_[xn_ind][xo_ind], nbf_);
timer_off("DDF pQq big K DGEMM");
			} else {
				X_Block_sparse('N', true, xn_iter, pruned_c, pruned_d, a, xn, u, nullptr, nullptr, eri);
				xn_ind = xn_iter;
timer_on("DDF pQq big K DGEMM");
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xo_ind], naux_*nocc, 1.0, xn, naux_*nocc, xo, naux_*nocc, 0.0, k + k_disps_[xn_ind][xo_ind], nbf_);
timer_off("DDF pQq big K DGEMM");
			}
		}
		xh = xn;
		xn = xo;//XO.get();
		xo = xh;//XN.get();
		xo_ind = xn_ind;
	}


	for (size_t kf_i = 0; kf_i < nbf_; kf_i++){
		for (size_t kf_j = 0; kf_j < kf_i; kf_j++){
			k[kf_i * nbf_ + kf_j] = k[ kf_j *nbf_ + kf_i];
		}
	}
	
	if (first_char != 'B') {
		X_Block_sparse('P', false, 0, pruned_c, pruned_d, a, nullptr, u, v, pruned_j, eri);
	}
}

void DirectDFJK::pQp_mn_sparse_set_mP() {
	mP_func_map_pQq_.clear();
	mP_shel_map_pQq_.clear();

	mP_func_map_pQq_.resize(p_shells_);
	mP_shel_map_pQq_.resize(p_shells_);
	
	prepare_p_blocks();

	printf("num_blocks_ is %zu\n", num_blocks_);

	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
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

	int nocc = C_left_ao_[0]->ncol();

	double* j = J_ao_[0]->pointer(0)[0];
	double* k = K_ao_[0]->pointer(0)[0];
	
	std::unique_ptr<double[]> A(new double[ biggest_shell_ * naux_ * nbf_ ]);
	std::unique_ptr<double[]> U(new double[biggest_shell_ * naux_*nocc]);
	std::unique_ptr<double[]> XN(new double[biggest_block_/nbf_*nocc]);
	std::unique_ptr<double[]> XO(new double[biggest_block_/nbf_*nocc]);
	std::unique_ptr<double[]> V(new double[naux_]);
	std::unique_ptr<double[]> P_V(new double[naux_]);
	std::unique_ptr<double[]> P_CM(new double[naux_*naux_]);


	double* a = A.get();
	double* u = U.get();
	double* xn = XN.get();
	double* xo = XO.get();
	double* xh;
	double* v = V.get();
	double* p_v = P_V.get();
	double* p_cm = P_CM.get();


	std::unique_ptr<double[]> P_C(new double[nbf_*nocc]);
	std::unique_ptr<double[]> P_D(new double[nbf_]);
	std::unique_ptr<double[]> P_J(new double[nbf_]);


	double* pruned_c = P_C.get();
	double* pruned_d = P_D.get();
	double* pruned_j = P_J.get();


	double Zero = 0.0;
	double* zero = &Zero;
	C_DCOPY( nbf_*nbf_, zero, 0, j, 1);
    C_DCOPY( nbf_*nbf_, zero, 0, k, 1);
	C_DCOPY(naux_, zero, 0, v, 1);

	size_t xn_ind;
	size_t xo_ind = 0;

	printf("finished initializing variables\n");

	printf("do_K_ %s\n", (do_K_ ? "true": "false") );

	X_Block_mn_sparse_set_mP( num_blocks_==1, do_K_, 0, pruned_c, pruned_d, a, xo, u, v, p_v, pruned_j, p_cm, eri);
	printf("done with x_block 0\n");
	C_DGEMM('N', 'T', Block_funcs_[0], Block_funcs_[0], naux_*nocc, 1.0, xo, naux_*nocc, xo, naux_*nocc, 1.0, k, nbf_);
	for ( size_t xo_iter = 0; xo_iter < num_blocks_; xo_iter++ ) {
		for (size_t xn_iter = 0; xn_iter < num_blocks_ - xo_iter ; xn_iter++ ) {
			if ( xo_iter == 0 && xn_iter < num_blocks_ - 1 ) {
				xn_ind = xn_iter;
				X_Block_mn_sparse_set_mP( false, do_K_, xn_ind, pruned_c, pruned_d, a, xn, u, v, p_v, pruned_j, p_cm, eri);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
			}
			if ( xo_iter == 0 && xn_iter == num_blocks_ - 1 ) {
				xn_ind = xn_iter;
				X_Block_mn_sparse_set_mP( true, do_K_, xn_ind, pruned_c, pruned_d, a, xn, u, v, p_v, pruned_j, p_cm, eri);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
			}
			if ( xo_iter == 1 ) {
				xn_ind = xn_iter;
				X_Block_mn_mP_sparse('P', do_K_, xn_ind, pruned_c, pruned_d, a, xn, u, v, p_v, pruned_j, p_cm, eri);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
			} else {
				xn_ind = xn_iter;
				X_Block_mn_mP_sparse('N', do_K_, xn_ind, pruned_c, pruned_d, a, xn, u, v, p_v, pruned_j, p_cm, eri);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
			}
		}
		xh = xn;
		xo = xn;
		xn = xh;
		xo_ind = xn_ind;
	}

	for (size_t kf_i = 0; kf_i < nbf_; kf_i++){
		for (size_t kf_j = 0; kf_j < kf_i; kf_j++){
			k[kf_i * nbf_ + kf_j] = k[ kf_j *nbf_ + kf_i];
		}
	}
	
	if (num_blocks_ != 1) {
		X_Block_sparse('P', false, 0, pruned_c, pruned_d, a, nullptr, u, v, pruned_j, eri);
	}

	nocc_last_ = nocc;

}

void DirectDFJK::pQp_mn_mP_sparse() {
	prepare_p_blocks();

	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
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

	int nocc = C_left_ao_[0]->ncol();

	double* j = J_ao_[0]->pointer(0)[0];
	double* k = K_ao_[0]->pointer(0)[0];
	
	std::unique_ptr<double[]> A(new double[ biggest_shell_ * naux_ * nbf_ ]);
	std::unique_ptr<double[]> U(new double[biggest_shell_ * naux_*nocc]);
	std::unique_ptr<double[]> XN(new double[biggest_block_/nbf_*nocc]);
	std::unique_ptr<double[]> XO(new double[biggest_block_/nbf_*nocc]);
	std::unique_ptr<double[]> V(new double[naux_]);
	std::unique_ptr<double[]> P_V(new double[naux_]);
	std::unique_ptr<double[]> P_CM(new double[naux_*naux_]);


	double* a = A.get();
	double* u = U.get();
	double* xn = XN.get();
	double* xo = XO.get();
	double* xh;
	double* v = V.get();
	double* p_v = P_V.get();
	double* p_cm = P_CM.get();


	std::unique_ptr<double[]> P_C(new double[nbf_*nocc]);
	std::unique_ptr<double[]> P_D(new double[nbf_]);
	std::unique_ptr<double[]> P_J(new double[nbf_]);


	double* pruned_c = P_C.get();
	double* pruned_d = P_D.get();
	double* pruned_j = P_J.get();


	double Zero = 0.0;
	double* zero = &Zero;
	C_DCOPY( nbf_*nbf_, zero, 0, j, 1);
    C_DCOPY( nbf_*nbf_, zero, 0, k, 1);
	C_DCOPY(naux_, zero, 0, v, 1);

	char first_char = ( num_blocks_ == 1 ? 'B' : 'V');

	size_t xn_ind;
	size_t xo_ind = 0;

	X_Block_mn_mP_sparse( first_char, do_K_, 0, pruned_c, pruned_d, a, xo, u, v, p_v, pruned_j, p_cm, eri);
	C_DGEMM('N', 'T', Block_funcs_[0], Block_funcs_[0], naux_*nocc, 1.0, xo, naux_*nocc, xo, naux_*nocc, 1.0, k, nbf_);
	for ( size_t xo_iter = 0; xo_iter < num_blocks_; xo_iter++ ) {
		for (size_t xn_iter = 0; xn_iter < num_blocks_ - xo_iter ; xn_iter++ ) {
			if ( xo_iter == 0 && xn_iter < num_blocks_ - 1 ) {
				xn_ind = xn_iter;
				X_Block_mn_mP_sparse('V' , do_K_, xn_ind, pruned_c, pruned_d, a, xn, u, v, p_v, pruned_j, p_cm, eri);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
			}
			if ( xo_iter == 0 && xn_iter == num_blocks_ - 1 ) {
				xn_ind = xn_iter;
				X_Block_mn_mP_sparse( 'B', do_K_, xn_ind, pruned_c, pruned_d, a, xn, u, v, p_v, pruned_j, p_cm, eri);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
				C_DGEMM('N', 'T', Block_funcs_[xn_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xn, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
			}
			if ( xo_iter == 1 ) {
				xn_ind = xn_iter;
				X_Block_mn_mP_sparse('P', do_K_, xn_ind, pruned_c, pruned_d, a, xn, u, v, p_v, pruned_j, p_cm, eri);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
			} else {
				xn_ind = xn_iter;
				X_Block_mn_mP_sparse('N', do_K_, xn_ind, pruned_c, pruned_d, a, xn, u, v, p_v, pruned_j, p_cm, eri);
				C_DGEMM('N', 'T', Block_funcs_[xo_ind], Block_funcs_[xn_ind], naux_*nocc, 1.0, xo, naux_*nocc, xn, naux_*nocc, 1.0, k, nbf_);
			}
		}
		xh = xn;
		xo = xn;
		xn = xh;
		xo_ind = xn_ind;
	}

	for (size_t kf_i = 0; kf_i < nbf_; kf_i++){
		for (size_t kf_j = 0; kf_j < kf_i; kf_j++){
			k[kf_i * nbf_ + kf_j] = k[ kf_j *nbf_ + kf_i];
		}
	}

	if (first_char != 'B') {
		X_Block_sparse('P', false, 0, pruned_c, pruned_d, a, nullptr, u, v, pruned_j, eri);
	}
}

// Function that produces tensors for the final contraction
//    of the exchange matrix build. The issue is that it may or may not
//    have to construct one of various terms for the coulomb matrix.
//    We handle this with a switch.
// coul_work \in { 'V', 'P', 'B', 'N' }
// 'V' means we compute a vector to be contracted against the coulomb Metric.
void DirectDFJK::X_Block( char coul_work, bool compute_k, size_t block, double* ao_block, double* x, double* u, double* coulomb_vector, std::vector<std::shared_ptr<TwoBodyAOInt>> eri){

	size_t nocc = C_left_ao_[0]->ncol();
	double* c = C_left_ao_[0]->pointer()[0];
	double* met_m_0_5 = get_metric_power(-0.5);
	double* j = J_ao_[0]->pointer()[0];
	double* d = D_ao_[0]->pointer()[0];

	std::unique_ptr<int[]> IPIV(new int[naux_]);
	int* ipiv = IPIV.get();
		
	switch (coul_work) {
		case 'N':
for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF AO_CONST");
	compute_dense_AO_block_p_pQq(shell_iter, ao_block, eri);
timer_off("DDF AO_CONST");
	if (compute_k) {
timer_on("DDF pQq small K DGEMM");
		C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, nbf_, 1.0, ao_block, nbf_, c, nocc, 0.0, u, nocc);
		for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
			C_DGEMM( 'N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_*nocc, nocc,  0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index() + func_it)* naux_*nocc, nocc);
		}
timer_off("DDF pQq small K DGEMM");
	}
}
			break;
		case 'V':
			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++){
				//compute ao blocks
timer_on("DDF AO_SPARSE");
				compute_dense_AO_block_p_pQq(shell_iter, ao_block , eri);
timer_off("DDF AO_SPARSE");
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
// Form V for Coulomb Matrix construction
					C_DGEMV( 'N', (int) naux_, (int) nbf_, 1.0, ao_block + func_it*naux_*nbf_, nbf_, d + (primary_->shell(shell_iter).function_index()+func_it)*nbf_, 1, 1.0, coulomb_vector, 1 );
				}
				if (compute_k) {
timer_on("DDF pQq small K DGEMM");
// Form U for Exchange Matrix construction
					C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, nbf_, 1.0, ao_block, nbf_, c, nocc, 0.0, u, nocc);
// Contract this u into the corresponding portion of x
					for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
						C_DGEMM( 'N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_*nocc, nocc,  0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index() + func_it)* naux_*nocc, nocc);
					}
timer_off("DDF pQq small K DGEMM");
				}
			}
			break;
		case 'P':
			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF AO_SPARSE");
				compute_dense_AO_block_p_pQq(shell_iter, ao_block, eri);
timer_off("DDF AO_SPARSE");
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
					C_DGEMV('T', (int) naux_, (int) nbf_, 1.0, ao_block + func_it * naux_ * nbf_ , nbf_, coulomb_vector, 1, 0.0, j + nbf_ * (primary_->shell(shell_iter).function_index() + func_it), 1 );
				}
				if (compute_k) {
timer_on("DDF pQq small K DGEMM");
					C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, nbf_, 1.0, ao_block, nbf_, c, nocc, 0.0, u, nocc);
					for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
						C_DGEMM( 'N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_*nocc, nocc,  0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index() + func_it)* naux_*nocc, nocc);
					}
timer_off("DDF pQq small K DGEMM");
				}
			}
			break;
		case 'B':
			double* metp = &CMPQ_LU_.front();
            int* pert = &PERMUTE_.front();
			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF AO_SPARSE");
            	compute_dense_AO_block_p_pQq(shell_iter, ao_block, eri);
timer_off("DDF AO_SPARSE");
timer_on("DDF pQq small K DGEMM");
				C_DGEMM( 'N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, nbf_, 1.0, ao_block, nbf_, c, nocc, 0.0, u, nocc );
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
					C_DGEMM( 'N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_*nocc, nocc,  0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index() + func_it)* naux_*nocc, nocc);
					C_DGEMV( 'N', naux_, nbf_, 1.0, ao_block +func_it*naux_*nbf_ , nbf_, d + ( primary_->shell(shell_iter).function_index()+ func_it)*nbf_, 1, 1.0, coulomb_vector, 1);
				}
timer_off("DDF pQq small K DGEMM");
			}

			C_DGETRS( 'N', naux_, 1, metp, naux_, pert, coulomb_vector, naux_);

			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
timer_on("DDF AO_SPARSE");
            	compute_dense_AO_block_p_pQq(shell_iter, ao_block, eri);
timer_off("DDF AO_SPARSE");
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
					C_DGEMV('T', naux_, nbf_, 1.0, ao_block + func_it * naux_ * nbf_, nbf_, coulomb_vector, 1, 0.0, j + ( primary_->shell(shell_iter).function_index() + func_it)*nbf_, 1 );
				}
			}
        break;
    }
}

void DirectDFJK::X_Block_sparse( char coul_work, bool compute_k, size_t block, double* pruned_c, double* pruned_d, double* ao_block, double* x, double* u, double* coulomb_vector, double* pruned_j, std::vector<std::shared_ptr<TwoBodyAOInt>> eri){

	size_t nocc = C_left_ao_[0]->ncol();
	double* c = C_left_ao_[0]->pointer()[0];
	double* met_m_0_5 = get_metric_power(-0.5);
	double* j = J_ao_[0]->pointer()[0];
	double* d = D_ao_[0]->pointer()[0];

	switch (coul_work) {
		case 'N':
			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
	compute_sparse_AO_block_p_pQq(shell_iter, ao_block, eri);
				if (compute_k) {
					prune_c( shell_iter, nocc, pruned_c, c );
					C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc);
					for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
						C_DGEMM( 'N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_*nocc, nocc,  0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index() + func_it)* naux_*nocc, nocc);
					}
				}
			}
			break;
		case 'V':
			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++){
				//compute ao blocks
	compute_sparse_AO_block_p_pQq(shell_iter, ao_block , eri);
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
					prune_d( shell_iter, pruned_d, d + nbf_ * (func_it + primary_->shell(shell_iter).function_index()) );
// Form V for Coulomb Matrix construction
					C_DGEMV( 'N', (int) naux_, (int) schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it*naux_*schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], pruned_d, 1, 1.0, coulomb_vector, 1 );
				}
				if (compute_k) {
					prune_c( shell_iter, nocc, pruned_c, c );
// Form U for Exchange Matrix construction
					C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc);
// Contract this u into the corresponding portion of x
					for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
						C_DGEMM( 'N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_*nocc, nocc,  0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index() + func_it)* naux_*nocc, nocc);
					}
				}
			}
			break;
		case 'P':
			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
	compute_sparse_AO_block_p_pQq(shell_iter, ao_block, eri);
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++) {
					C_DGEMV('T', (int) naux_, (int) schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it * naux_ * schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], coulomb_vector, 1, 0.0, pruned_j, 1 );
					unprune_J( shell_iter, j + nbf_ * (primary_->shell(shell_iter).function_index() + func_it), pruned_j);
				}
				if (compute_k) {
					prune_c( shell_iter, nocc, pruned_c, c );
					C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc);
					for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
						C_DGEMM( 'N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_*nocc, nocc,  0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index() + func_it)* naux_*nocc, nocc);
					}
				}
			}
			break;
		case 'B':
			double* metp = &CMPQ_LU_.front();
            int* pert = &PERMUTE_.front();
			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
    compute_sparse_AO_block_p_pQq(shell_iter, ao_block, eri);
				if (compute_k) {
					prune_c( shell_iter, nocc, pruned_c, c );
					C_DGEMM( 'N', 'N', primary_->shell(shell_iter).nfunction()*naux_, nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc );
					for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
						C_DGEMM( 'N', 'N', naux_, nocc, naux_, 1.0, met_m_0_5, naux_, u + func_it * naux_*nocc, nocc,  0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(Shell_starts_[block]).function_index() + func_it)* naux_*nocc, nocc);
					}
				}
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
					prune_d( shell_iter, pruned_d, d + nbf_ * (func_it + primary_->shell(shell_iter).function_index()) );
					C_DGEMV( 'N', naux_, schwarz_dense_funcs_[shell_iter], 1.0, ao_block +func_it*naux_*schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], pruned_d, 1, 1.0, coulomb_vector, 1);
				}
			}


			C_DGETRS( 'N', naux_, 1, metp, naux_, pert, coulomb_vector, naux_);
			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
            	compute_sparse_AO_block_p_pQq(shell_iter, ao_block, eri);
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
					C_DGEMV( 'T', (int) naux_, (int) schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it * naux_ * schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], coulomb_vector, 1, 0.0, pruned_j, 1 );
					unprune_J( shell_iter, j + nbf_ * (primary_->shell(shell_iter).function_index() + func_it), pruned_j);
				}
			}
        break;
    }
}

void DirectDFJK::X_Block_mn_sparse_set_mP(bool with_contraction, bool compute_k, size_t block, double* pruned_c, double* pruned_d, double* ao_block, double* x, double* u, double* coulomb_vector, double* pruned_coulomb_vector, double* pruned_j, double* pruned_cm, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ){

    size_t nocc = C_left_ao_[0]->ncol();
    double* c = C_left_ao_[0]->pointer()[0];
    double* met_m_0_5 = get_metric_power(-0.5);
    double* j = J_ao_[0]->pointer()[0];
    double* d = D_ao_[0]->pointer()[0];

//!NOT!NOT!NOT!NOT!NOT!NOT!NOT! with_contraction
    if (!with_contraction){
		printf("in X_Block set_mP without contraction\n");
        for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
            compute_AO_block_p_pQq_mn_sparse_set_mP( shell_iter, ao_block, eri );
            if (compute_k) {
                prune_c( shell_iter, nocc, pruned_c, c );
     			prune_cmpq( shell_iter, met_m_0_5, pruned_cm);
				C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * mP_func_map_pQq_[shell_iter].size(), nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc );
				for ( size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
					C_DGEMM('N', 'N', naux_, nocc, mP_func_map_pQq_[shell_iter].size(), 1.0, pruned_cm, mP_func_map_pQq_[shell_iter].size(), u + func_it * nocc * mP_func_map_pQq_[shell_iter].size(), nocc, 0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(shell_iter).function_index() + func_it)*nocc*naux_, nocc );
				}
            }
			for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
				prune_d( shell_iter, pruned_d, d + nbf_ * (func_it + primary_->shell(shell_iter).function_index()) );
				C_DGEMV( 'N', mP_func_map_pQq_[shell_iter].size(), schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it * mP_func_map_pQq_[shell_iter].size() * schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], pruned_d, 1, 1.0, pruned_coulomb_vector, 1); 
				unprune_V( shell_iter, coulomb_vector , pruned_coulomb_vector);
			}
		}
	} else {
		printf("in X_Block set_mP with contraction\n");
        double* metp = &CMPQ_LU_.front();
        int* pert = &PERMUTE_.front();
		printf("got arrays\n");
        for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
			printf("about to get aos\n");
            compute_AO_block_p_pQq_mn_mP_sparse( shell_iter, ao_block, eri );
			printf("got aos\n");
			if (compute_k) {
				printf("compute_k true\n");
            	prune_c( shell_iter, nocc, pruned_c, c );
				prune_cmpq( shell_iter, met_m_0_5, pruned_cm);
				C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * mP_func_map_pQq_[shell_iter].size(), nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc );
				for ( size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
					C_DGEMM('N', 'N', naux_, nocc, mP_func_map_pQq_[shell_iter].size(), 1.0, pruned_cm, mP_func_map_pQq_[shell_iter].size(), u + func_it * nocc * mP_func_map_pQq_[shell_iter].size(), nocc, 0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(shell_iter).function_index() + func_it)*nocc*naux_, nocc );
				}
			}
			for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
				prune_d( shell_iter, pruned_d, d + nbf_ * (func_it + primary_->shell(shell_iter).function_index()) );
				C_DGEMV( 'N', mP_func_map_pQq_[shell_iter].size(), schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it * mP_func_map_pQq_[shell_iter].size() * schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], pruned_d, 1, 1.0, pruned_coulomb_vector, 1); 
				unprune_V( shell_iter, coulomb_vector , pruned_coulomb_vector);
			}
        }
		printf("about to solve JPHI=V\n");
        C_DGETRS( 'N', naux_, 1, metp, naux_, pert, coulomb_vector, naux_);
		printf("solved JPHI=V\n");
        for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
            compute_AO_block_p_pQq_mn_mP_sparse( shell_iter, ao_block, eri );
			printf("about to prune phi\n");
			prune_phi( shell_iter, coulomb_vector, pruned_coulomb_vector  );
			printf("pruned phi\n");
			for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
				printf("about to contract into J\n");
				C_DGEMV( 'T', (int) mP_shel_map_pQq_[shell_iter].size(), (int) schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it * mP_shel_map_pQq_[shell_iter].size() * schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], coulomb_vector, 1, 0.0, pruned_j, 1 );
				printf("contracted into pruned_j\n");
				printf("unpruning j\n");
				unprune_J( shell_iter, j + nbf_ * (primary_->shell(shell_iter).function_index() + func_it), pruned_j);
				printf("j unpruned\n");
			}
        }
    }

}


void DirectDFJK::X_Block_mn_mP_sparse(char coul_work, bool compute_k, size_t block, double* pruned_c, double* pruned_d, double* ao_block, double* x, double* u, double* coulomb_vector, double* pruned_coulomb_vector, double* pruned_j, double* pruned_cm, std::vector<std::shared_ptr<TwoBodyAOInt>> eri ){

    size_t nocc = C_left_ao_[0]->ncol();
    double* c = C_left_ao_[0]->pointer()[0];
    double* met_m_0_5 = get_metric_power(-0.5);
    double* j = J_ao_[0]->pointer()[0];
    double* d = D_ao_[0]->pointer()[0];


    switch (coul_work) {
        case 'N':
        	for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
            	compute_AO_block_p_pQq_mn_mP_sparse( shell_iter, ao_block, eri );
        	    if (compute_k) {
            	    prune_c( shell_iter, nocc, pruned_c, c );
     				prune_cmpq( shell_iter, met_m_0_5, pruned_cm);
					C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * mP_func_map_pQq_[shell_iter].size(), nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc );
					for ( size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
						C_DGEMM('N', 'N', naux_, nocc, mP_func_map_pQq_[shell_iter].size(), 1.0, pruned_cm, mP_func_map_pQq_[shell_iter].size(), u + func_it * nocc * mP_func_map_pQq_[shell_iter].size(), nocc, 0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(shell_iter).function_index() + func_it)*nocc*naux_, nocc );
					}
            	}
			}
            break;
        case 'V':
        	for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
            	compute_AO_block_p_pQq_mn_mP_sparse( shell_iter, ao_block, eri );
            	if (compute_k) {
                	prune_c( shell_iter, nocc, pruned_c, c );
     				prune_cmpq( shell_iter, met_m_0_5, pruned_cm);
					C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * mP_func_map_pQq_[shell_iter].size(), nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc );
					for ( size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
						C_DGEMM('N', 'N', naux_, nocc, mP_func_map_pQq_[shell_iter].size(), 1.0, pruned_cm, mP_func_map_pQq_[shell_iter].size(), u + func_it * nocc * mP_func_map_pQq_[shell_iter].size(), nocc, 0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(shell_iter).function_index() + func_it)*nocc*naux_, nocc );
					}
            	}
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
					prune_d( shell_iter, pruned_d, d + nbf_ * (func_it + primary_->shell(shell_iter).function_index()) );
					C_DGEMV( 'N', mP_func_map_pQq_[shell_iter].size(), schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it * mP_func_map_pQq_[shell_iter].size() * schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], pruned_d, 1, 1.0, pruned_coulomb_vector, 1); 
					unprune_V( shell_iter, coulomb_vector , pruned_coulomb_vector);
				}
			}
            break;
        case 'P':
	        for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
  		        compute_AO_block_p_pQq_mn_mP_sparse( shell_iter, ao_block, eri );
        	    if (compute_k) {
            	    prune_c( shell_iter, nocc, pruned_c, c );
     				prune_cmpq( shell_iter, met_m_0_5, pruned_cm);
					C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * mP_func_map_pQq_[shell_iter].size(), nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc );
					for ( size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
						C_DGEMM('N', 'N', naux_, nocc, mP_func_map_pQq_[shell_iter].size(), 1.0, pruned_cm, mP_func_map_pQq_[shell_iter].size(), u + func_it * nocc * mP_func_map_pQq_[shell_iter].size(), nocc, 0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(shell_iter).function_index() + func_it)*nocc*naux_, nocc );
					}
            	}
			}	
			for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
				prune_phi( shell_iter, coulomb_vector, pruned_coulomb_vector  );
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
					C_DGEMV( 'T', (int) mP_shel_map_pQq_[shell_iter][func_it], (int) schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it * mP_shel_map_pQq_[shell_iter][func_it] * schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], coulomb_vector, 1, 0.0, pruned_j, 1 );
					unprune_J( shell_iter, j + nbf_ * (primary_->shell(shell_iter).function_index() + func_it), pruned_j);
				}
			}
            break;
        case 'B':
        	double* metp = &CMPQ_LU_.front();
        	int* pert = &PERMUTE_.front();
        	for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
            	compute_AO_block_p_pQq_mn_mP_sparse( shell_iter, ao_block, eri );
				if (compute_k) {
            		prune_c( shell_iter, nocc, pruned_c, c );
					prune_cmpq( shell_iter, met_m_0_5, pruned_cm);
					C_DGEMM('N', 'N', primary_->shell(shell_iter).nfunction() * mP_func_map_pQq_[shell_iter].size(), nocc, schwarz_dense_funcs_[shell_iter], 1.0, ao_block, schwarz_dense_funcs_[shell_iter], pruned_c, nocc, 0.0, u, nocc );
					for ( size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
						printf("%zu %zu %zu\n", Q_shells_, naux_, mP_func_map_pQq_[shell_iter].size());
						C_DGEMM('N', 'N', naux_, nocc, mP_func_map_pQq_[shell_iter].size(), 1.0, pruned_cm, mP_func_map_pQq_[shell_iter].size(), u + func_it * nocc * mP_func_map_pQq_[shell_iter].size(), nocc, 0.0, x + (primary_->shell(shell_iter).function_index() - primary_->shell(shell_iter).function_index() + func_it)*nocc*naux_, nocc );
					}
				}
				for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++ ){
					prune_d( shell_iter, pruned_d, d + nbf_ * (func_it + primary_->shell(shell_iter).function_index()) );
					C_DGEMV( 'N', mP_func_map_pQq_[shell_iter].size(), schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it * mP_func_map_pQq_[shell_iter].size() * schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], pruned_d, 1, 1.0, pruned_coulomb_vector, 1); 
					unprune_V( shell_iter, coulomb_vector , pruned_coulomb_vector);
				}
        	}
        	C_DGETRS( 'N', naux_, 1, metp, naux_, pert, coulomb_vector, naux_);
        	for (size_t shell_iter = Shell_starts_[block]; shell_iter <= Shell_stops_[block]; shell_iter++) {
            	compute_AO_block_p_pQq_mn_mP_sparse( shell_iter, ao_block, eri );
				prune_phi( shell_iter, coulomb_vector, pruned_coulomb_vector  );
					for (size_t func_it = 0; func_it < primary_->shell(shell_iter).nfunction(); func_it++){
						C_DGEMV( 'T', (int) mP_shel_map_pQq_[shell_iter][func_it], (int) schwarz_dense_funcs_[shell_iter], 1.0, ao_block + func_it * mP_shel_map_pQq_[shell_iter][func_it] * schwarz_dense_funcs_[shell_iter], schwarz_dense_funcs_[shell_iter], coulomb_vector, 1, 0.0, pruned_j, 1 );
						unprune_J( shell_iter, j + nbf_ * (primary_->shell(shell_iter).function_index() + func_it), pruned_j);
					}
        	}
            break;
    }

}

void DirectDFJK::prune_c( size_t &mu, size_t nocc, double* pruned_c, double* raw_c ) {
	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
#endif
	size_t pru_add;
	size_t raw_add;

#pragma omp parallel for schedule(dynamic) num_threads(procs)
	for (size_t shell_it = 0U; shell_it < schwarz_shell_mask_pQq_[mu].size(); shell_it++){
		for (size_t in_shell = 0U; in_shell < primary_->shell(schwarz_shell_mask_pQq_[mu][shell_it]).nfunction(); in_shell++ ) {
			for (size_t mol_it = 0U; mol_it < nocc; mol_it++) {

				pruned_c[ ( schwarz_func_starts_pQq_[mu][shell_it] + in_shell ) * nocc+ mol_it ] 
				= 
				raw_c[(primary_->shell( schwarz_shell_mask_pQq_[mu][shell_it]).function_index() + in_shell)*nocc + mol_it ];
			}
		}
	}

}

void DirectDFJK::prune_d( size_t &mu, double* pruned_d, double* raw_d){
	int procs = 1;
#ifdef _OPENMP
	procs = omp_nthread_;
#endif

	for (size_t shell_it = 0U; shell_it < schwarz_shell_mask_pQq_[mu].size(); shell_it++) {
		for (size_t in_shell = 0U; in_shell < primary_->shell(schwarz_shell_mask_pQq_[mu][shell_it]).nfunction(); in_shell++) {
			pruned_d[ schwarz_func_starts_pQq_[mu][shell_it] + in_shell ]
			=
			raw_d[primary_->shell(schwarz_shell_mask_pQq_[mu][shell_it]).function_index() + in_shell];
								  
		}
	}
	
}

// You should only pass in one row of J in at a time. We will trust
//    the wrapping function to take care of this.
// mu is a shell index. we need it to get at sparsity information.
void DirectDFJK::unprune_J( size_t &mu, double* raw_j, double* pruned_j ) {
	size_t address = 0;
	for (size_t shell_iter = 0; shell_iter < schwarz_shell_mask_pQq_[mu].size(); shell_iter++) {
		for (size_t in_shell = 0; in_shell < primary_->shell( schwarz_shell_mask_pQq_[mu][shell_iter] ).nfunction(); in_shell++ ) {
			raw_j[primary_->shell( schwarz_shell_mask_pQq_[mu][shell_iter]).function_index() + in_shell] = 
			pruned_j[schwarz_func_starts_pQq_[mu][shell_iter] + in_shell];
		}
	}
}

void DirectDFJK::prune_cmpq(size_t big_Mu, double* raw_CMPQ, double* pruned_CMPQ) {

#pragma omp parallel for
    for ( size_t row_iter = 0; row_iter < naux_; row_iter++ ) {
        for ( size_t col_iter = 0; col_iter < mP_func_map_pQq_[big_Mu].size(); col_iter++ ) {
            pruned_CMPQ[ row_iter * mP_func_map_pQq_[big_Mu].size() + col_iter] = 
			raw_CMPQ[ row_iter * naux_ + mP_func_map_pQq_[big_Mu][col_iter]] ;
        }
    }
}

void DirectDFJK::unprune_V( size_t big_Mu, double* raw_v, double* pruned_v){
	for (size_t i = 0; i < mP_func_map_pQq_[big_Mu].size(); i++ ) {
		raw_v[mP_func_map_pQq_[big_Mu][i]] += pruned_v[i];
	}
}

void DirectDFJK::prune_phi( size_t big_Mu, double* raw_phi, double* pruned_phi) {
	for ( size_t i = 0; i < mP_func_map_pQq_[big_Mu].size(); i++ ) {
		pruned_phi[i] = raw_phi[ mP_func_map_pQq_[big_Mu][i] ];
	}
}

} //namespace psi
