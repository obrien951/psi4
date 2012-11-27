// p_so(pitzer) = p_symblk + PitzerOffset[h]; where h=mosym[p_symblk]
// p_symblk = pitzer2symblk[p_so(pitzer)];

/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string> 
#include <iomanip>
#include <vector> 


/** Required PSI3 includes */ 
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>

/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "occwave.h"
#include "defines.h"
#include "arrays.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace occwave{

void OCCWave::get_moinfo()
{      
  //fprintf(outfile,"\n get_moinfo is starting... \n"); fflush(outfile);
/********************************************************************************************/
/************************** MO info *********************************************************/
/********************************************************************************************/
	// Read in mo info
        nso_     = reference_wavefunction_->nso();
        nirrep_ = reference_wavefunction_->nirrep();
        nmo_     = reference_wavefunction_->nmo();
        nmopi_    = reference_wavefunction_->nmopi();
        nsopi_    = reference_wavefunction_->nsopi();
        doccpi_  = reference_wavefunction_->doccpi();
        soccpi_  = reference_wavefunction_->soccpi();
        frzcpi_  = reference_wavefunction_->frzcpi();
        frzvpi_  = reference_wavefunction_->frzvpi();

        // get nfrzc and nfrzv
        nfrzc = 0;
        nfrzv = 0;
        for(int h=0; h<nirrep_; h++) {
	  nfrzc += frzcpi_[h];
	  nfrzv += frzvpi_[h];
	}
	
	// form occpi and virtpi
	occpiA = init_int_array(nirrep_);
	virtpiA = init_int_array(nirrep_);
	memset(occpiA,0, sizeof(int)*nirrep_);
	memset(virtpiA,0, sizeof(int)*nirrep_);
	for(int h=0; h<nirrep_; h++) {
	  virtpiA[h] = nmopi_[h] - doccpi_[h];
	  occpiA[h] = doccpi_[h];
	}
	
	//active occ and virt
	adoccpi = init_int_array(nirrep_);
	aoccpiA = init_int_array(nirrep_);
	avirtpiA = init_int_array(nirrep_);
	memset(adoccpi,0, sizeof(int)*nirrep_);
	memset(aoccpiA,0, sizeof(int)*nirrep_);
	memset(avirtpiA,0, sizeof(int)*nirrep_);
	for(int h=0; h<nirrep_; h++) {
	  adoccpi[h] = doccpi_[h] - frzcpi_[h];
	  avirtpiA[h] = virtpiA[h] - frzvpi_[h];
	  aoccpiA[h] = doccpi_[h] - frzcpi_[h];
	}
	
	// Read in nuclear repulsion energy
	Enuc = Process::environment.molecule()->nuclear_repulsion_energy();
	
	// Read SCF energy
        Escf=reference_wavefunction_->reference_energy();
	Eref=Escf;
	Eelec=Escf-Enuc;
	
	/* Build mosym arrays */
	mosym = new int [nmo_];
	memset(mosym,0,sizeof(int)*nmo_);
	for(int h=0, q=0; h < nirrep_; h++){
	  for(int p=0; p < nmopi_[h]; p++){
	    mosym[q++] = h;
	  }
	}
	
	/* Build sosym arrays */
	sosym = new int [nso_];
	memset(sosym,0,sizeof(int)*nmo_);
	for(int h=0, q=0; h < nirrep_; h++){
	  for(int p=0; p < nsopi_[h]; p++){
	    sosym[q++] = h;
	  }
	}

	// find nooA
	nooA=0;
	for(int h=0; h < nirrep_; h++){
	  for(int p=0; p < doccpi_[h]; p++){
	    nooA++;
	  }
	}

	// PitzerOffset 
	PitzerOffset = new int[nirrep_];
	memset(PitzerOffset,0,sizeof(int)*nirrep_);
	for(int h=1; h < nirrep_; h++){
	  PitzerOffset[h] = PitzerOffset[h-1] + nmopi_[h-1];
	}
	
	nvoA=nmo_-nooA;   	// Number of virtual orbitals
	nacooA=nooA-nfrzc; 	// Number of active occupied orbitals
	nacso=nmo_-nfrzc-nfrzv; 	// Number of active  orbitals
	nacvoA=nvoA-nfrzv; 	// Number of active virtual orbitals
	npop=nmo_-nfrzv;         // Number of populated orbitals

	ntri_so = 0.5*nso_*(nso_+1);
        ntri = 0.5*nmo_*(nmo_+1);
	dimtei = 0.5*ntri*(ntri+1);

/********************************************************************************************/
/************************** pitzer2symblk ***************************************************/
/********************************************************************************************/
      pitzer2symirrep = new int[nmo_];
      pitzer2symblk = new int[nmo_];
      occ2symblkA = new int[nooA];
      virt2symblkA = new int[nvoA];
      memset(pitzer2symirrep,0,sizeof(int)*nmo_);
      memset(pitzer2symblk,0,sizeof(int)*nmo_);
      memset(occ2symblkA,0,sizeof(int)*nooA);
      memset(virt2symblkA,0,sizeof(int)*nvoA);

      // pitzer2symblk
      int ij,myoffset;
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<nmopi_[h]; ++i) {
            pitzer2symirrep[ij] = h;
            pitzer2symblk[ij] = ij-myoffset;
            ij++;
        }
        myoffset += nmopi_[h];
      }
      
      // occ2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<occpiA[h]; ++i) {
            occ2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += occpiA[h];
      }
      
      
      // vir2symblkA
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<virtpiA[h]; ++i) {
            virt2symblkA[ij] = ij-myoffset;
            ij++;
        }
        myoffset += virtpiA[h];
      }

      
      // print
      if (print_ > 1) {
          for(int p = 0; p < nmo_; p++) {
              fprintf(outfile," p, pitzer2symblk[p]: %2d %2d \n", p, pitzer2symblk[p]);
              fflush(outfile);
          }
          fprintf(outfile,"\n");
          fflush(outfile);
      }
    
/********************************************************************************************/
/************************** qt2pitzer *******************************************************/
/********************************************************************************************/
      qt2pitzerA = new int[nmo_];
      pitzer2qtA = new int[nmo_];
      memset(qt2pitzerA,0,sizeof(int)*nmo_);
      memset(pitzer2qtA,0,sizeof(int)*nmo_);
 
if (reference_ == "RESTRICTED") {
      reorder_qt(doccpi_, soccpi_, frzcpi_, frzvpi_, pitzer2qtA, nmopi_, nirrep_);
      for(int p = 0; p < nmo_; p++) {
	  int pa = pitzer2qtA[p];
	  qt2pitzerA[pa] = p;
      }

      // print
      if (print_ > 1) {
      for(int p = 0; p < nmo_; p++) {
          fprintf(outfile," p, pitzer2qtA[p]: %2d %2d \n", p, pitzer2qtA[p]);
          fflush(outfile);
      }
         fprintf(outfile,"\n");
         fflush(outfile);

      for(int p = 0; p < nmo_; p++) {
          fprintf(outfile," p, qt2pitzerA[p]: %2d %2d \n", p, qt2pitzerA[p]);
          fflush(outfile);
      }
         fprintf(outfile,"\n");
         fflush(outfile);
      }
}// end if (reference_ == "RESTRICTED") 
 
/********************************************************************************************/
/************************** occ_off & vir_off ***********************************************/
/********************************************************************************************/ 
    // occ_qt = occ_sym_block + occ_off =>convert occ sym block index to occ qt index
    // general_qt = occ_sym_block + occ_off => convert occ sym block index to general qt index 
    // vir_qt = vir_sym_block + vir_off 
    // general_qt = vir_sym_block + vir_off + nocc
    // gen_qt = occ_qt => for occupieds
    // gen_qt = vir_qt + nocc => for virtuals
    occ_offA = new int[nirrep_];
    vir_offA = new int[nirrep_];
    memset(occ_offA, 0, sizeof(int)*nirrep_);
    memset(vir_offA, 0, sizeof(int)*nirrep_);
    int ocountA = occpiA[0]; 
    int vcountA = virtpiA[0];
    for(int h=1; h < nirrep_; h++) {
      occ_offA[h] = ocountA;
      ocountA += occpiA[h];
      vir_offA[h] = vcountA;
      vcountA += virtpiA[h];
    }

      // print
      if (print_ > 1) {
          for(int h = 0; h < nirrep_; h++) {
              fprintf(outfile," h, occ_offA[h]: %2d %2d \n", h, occ_offA[h]);
              fflush(outfile);
          }
          fprintf(outfile,"\n");
          fflush(outfile);

          for(int h = 0; h < nirrep_; h++) {
              fprintf(outfile," h, vir_offA[h]: %2d %2d \n", h, vir_offA[h]);
              fflush(outfile);
          }
          fprintf(outfile,"\n");
          fflush(outfile);
      }
      
/********************************************************************************************/
/************************** pairs per irrep *************************************************/
/********************************************************************************************/
        oo_pairpiAA = new int[nirrep_];
        ov_pairpiAA = new int[nirrep_];
        vv_pairpiAA = new int[nirrep_];
        memset(oo_pairpiAA,0,sizeof(int)*nirrep_);
        memset(ov_pairpiAA,0,sizeof(int)*nirrep_);
        memset(vv_pairpiAA,0,sizeof(int)*nirrep_);
        for (int h1 = 0; h1 < nirrep_; h1++) {
            for (int h2 = 0; h2 < nirrep_; h2++) {
                 int h = h1^h2;
                 oo_pairpiAA[h] += occpiA[h1] * occpiA[h2];
                 ov_pairpiAA[h] += occpiA[h1] * virtpiA[h2];
                 vv_pairpiAA[h] += virtpiA[h1] * virtpiA[h2];
            }
        }

        if (print_ > 1) {
         for(int h=0; h < nirrep_; h++) {
            fprintf(outfile," h, oo_pairpiAA[h]: %2d %2d \n", h, oo_pairpiAA[h]);
            fflush(outfile);
         }
         fprintf(outfile,"\n");
         fflush(outfile);

         for(int h=0; h < nirrep_; h++) {
            fprintf(outfile," h, ov_pairpiAA[h]: %2d %2d \n", h, ov_pairpiAA[h]);
            fflush(outfile);
         }
         fprintf(outfile,"\n");
         fflush(outfile);

         for(int h=0; h < nirrep_; h++) {
            fprintf(outfile," h, vv_pairpiAA[h]: %2d %2d \n", h, vv_pairpiAA[h]);
            fflush(outfile);
         }
         fprintf(outfile,"\n");
         fflush(outfile);
        }

/********************************************************************************************/
/************************** pair indices ****************************************************/
/********************************************************************************************/
        int *itemppi = new int[nirrep_];

        // OO-pair
        oo_pairidxAA = new Array3i("oo_pairidxAA", nirrep_, nooA, nooA);
        oo_pairidxAA->zero();
        memset(itemppi,0,sizeof(int)*nirrep_);
        for (int h1 = 0; h1 < nirrep_; h1++) {
            for (int h2 = 0; h2 < nirrep_; h2++) {
                 int h = h1^h2;
                 for (int i = 0; i < occpiA[h1]; i++) {
                      int I = i + occ_offA[h1];
                      for (int j = 0; j < occpiA[h2]; j++) {
                           int J = j + occ_offA[h2];
                           oo_pairidxAA->set(h, I, J, itemppi[h]);
                           itemppi[h]++;
                      }
                 }
            }
        }
        if (print_ > 1) oo_pairidxAA->print();

        // VV pair
        vv_pairidxAA = new Array3i("vv_pairidxAA", nirrep_, nvoA, nvoA);
        vv_pairidxAA->zero();
        memset(itemppi,0,sizeof(int)*nirrep_);
        for (int h1 = 0; h1 < nirrep_; h1++) {
            for (int h2 = 0; h2 < nirrep_; h2++) {
                 int h = h1^h2;
                 int pcount = 0;
                 for (int a = 0; a < virtpiA[h1]; a++) {
                      int A = a + vir_offA[h1];
                      for (int b = 0; b < virtpiA[h2]; b++) {
                           int B = b + vir_offA[h2];
                           vv_pairidxAA->set(h, A, B, itemppi[h]);
                           itemppi[h]++;
                      }
                 }
            }
        }
        if (print_ > 1) vv_pairidxAA->print();
        delete [] itemppi;

/********************************************************************************************/
/************************** Read orbital coefficients ***************************************/
/********************************************************************************************/
        // read orbital coefficients from chkpt
	Ca_ = SharedMatrix(reference_wavefunction_->Ca());
	Ca_ref = boost::shared_ptr<Matrix>(new Matrix("Ref alpha MO coefficients", nirrep_, nsopi_, nmopi_));
	
	// read orbital coefficients from external files
	if (read_mo_coeff == "TRUE"){
	  fprintf(outfile,"\n\tReading MO coefficients in pitzer order from the external file CmoA.psi...\n");  
	  fflush(outfile);
	  double **C_pitzerA = block_matrix(nso_,nmo_);
	  memset(C_pitzerA[0], 0, sizeof(double)*nso_*nmo_);
	
	  // read binary data
	  ifstream InFile1;
	  InFile1.open("CmoA.psi", ios::in | ios::binary);
	  InFile1.read( (char*)C_pitzerA[0], sizeof(double)*nso_*nmo_);
	  InFile1.close();
	  
	  //set C_scf
	  Ca_->set(C_pitzerA);
	  free_block(C_pitzerA);
        }
        
        // Build Reference MOs
        Ca_ref->copy(Ca_);
	if(print_ > 2) Ca_->print();

/********************************************************************************************/
/************************** Create all required matrice *************************************/
/********************************************************************************************/
        // Build Hso
	Hso = boost::shared_ptr<Matrix>(new Matrix("SO-basis One-electron Ints", nirrep_, nsopi_, nsopi_));
	Tso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Kinetic Energy Ints", nirrep_, nsopi_, nsopi_));
	Vso = boost::shared_ptr<Matrix>(new Matrix("SO-basis Potential Energy Ints", nirrep_, nsopi_, nsopi_));
	Hso->zero();
	Tso->zero();
	Vso->zero();
	
	// Read SO-basis one-electron integrals
	double *so_ints = init_array(ntri_so);
        IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_T, so_ints, ntri_so, 0, 0, outfile);
        Tso->set(so_ints);
        IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_V, so_ints, ntri_so, 0, 0, outfile);
        Vso->set(so_ints);
        free(so_ints);
	Hso->copy(Tso); 
	Hso->add(Vso);
  
//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
if (reference_ == "UNRESTRICTED") {

/********************************************************************************************/
/************************** MO info *********************************************************/
/********************************************************************************************/
	// form occpi and virtpi
	occpiB = init_int_array(nirrep_);
	virtpiB = init_int_array(nirrep_);
	memset(occpiB,0, sizeof(int)*nirrep_);
	memset(virtpiB,0, sizeof(int)*nirrep_);
	for(int h=0; h<nirrep_; h++) {
	  virtpiB[h] = nmopi_[h] - doccpi_[h];
	  occpiB[h] = doccpi_[h];
	}
	
	//active occ and virt
	aoccpiB = init_int_array(nirrep_);
	avirtpiB = init_int_array(nirrep_);
	memset(aoccpiB,0, sizeof(int)*nirrep_);
	memset(avirtpiB,0, sizeof(int)*nirrep_);
	for(int h=0; h<nirrep_; h++) {
	  avirtpiB[h] = virtpiB[h] - frzvpi_[h];
	  aoccpiB[h] = doccpi_[h] - frzcpi_[h];
	}

	// find nooB
	nooB=0;
	for(int h=0; h < nirrep_; h++){
	  for(int p=0; p < doccpi_[h]; p++){
	    nooB++;
	  }
	}
	
	nooA=nooB;
	for(int h=0; h < nirrep_; h++){
	  for(int p=0; p < soccpi_[h]; p++){
	    nooA++;
	  }
	}
	
	nvoB=nmo_-nooB;   	// Number of virtual orbitals
	nacooB=nooB-nfrzc; 	// Number of active occupied orbitals
	nacvoB=nvoB-nfrzv; 	// Number of active virtual orbitals

/********************************************************************************************/
/************************** pitzer2symblk ***************************************************/
/********************************************************************************************/
      occ2symblkB = new int[nooB];
      virt2symblkB = new int[nvoB];
      memset(occ2symblkB,0,sizeof(int)*nooB);
      memset(virt2symblkB,0,sizeof(int)*nvoB);

      // occ2symblkB
      int ij,myoffset;
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<occpiB[h]; ++i) {
            occ2symblkB[ij] = ij-myoffset;
            ij++;
        }
        myoffset += occpiB[h];
      }
      
      // vir2symblkB
      ij = 0;
      myoffset = 0;
      for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<virtpiB[h]; ++i) {
            virt2symblkB[ij] = ij-myoffset;
            ij++;
        }
        myoffset += virtpiB[h];
      }

/********************************************************************************************/
/************************** qt2pitzer *******************************************************/
/********************************************************************************************/
      qt2pitzerB = new int[nmo_];
      pitzer2qtB = new int[nmo_];
      memset(qt2pitzerB,0,sizeof(int)*nmo_);
      memset(pitzer2qtB,0,sizeof(int)*nmo_);
      reorder_qt_uhf(doccpi_, soccpi_, frzcpi_, frzvpi_, pitzer2qtA, pitzer2qtB, nmopi_, nirrep_);
      for(int p = 0; p < nmo_; p++) {
	  int pa = pitzer2qtA[p];
	  int pb = pitzer2qtB[p];
	  qt2pitzerA[pa] = p;
	  qt2pitzerB[pb] = p;
      }

       // print
      if (print_ > 1) {
      for(int p = 0; p < nmo_; p++) {
          fprintf(outfile," p, pitzer2qtA[p]: %2d %2d \n", p, pitzer2qtA[p]);
          fflush(outfile);
      }
         fprintf(outfile,"\n");
         fflush(outfile);

      for(int p = 0; p < nmo_; p++) {
          fprintf(outfile," p, pitzer2qtB[p]: %2d %2d \n", p, pitzer2qtB[p]);
          fflush(outfile);
      }
         fprintf(outfile,"\n");
         fflush(outfile);

      for(int p = 0; p < nmo_; p++) {
          fprintf(outfile," p, qt2pitzerA[p]: %2d %2d \n", p, qt2pitzerA[p]);
          fflush(outfile);
      }
         fprintf(outfile,"\n");
         fflush(outfile);

      for(int p = 0; p < nmo_; p++) {
          fprintf(outfile," p, qt2pitzerB[p]: %2d %2d \n", p, qt2pitzerB[p]);
          fflush(outfile);
      }
         fprintf(outfile,"\n");
         fflush(outfile);
      }// end if   

/********************************************************************************************/
/************************** occ_off & vir_off ***********************************************/
/********************************************************************************************/ 
    occ_offB = new int[nirrep_];
    vir_offB = new int[nirrep_];
    memset(occ_offB, 0, sizeof(int)*nirrep_);
    memset(vir_offB, 0, sizeof(int)*nirrep_);
    int ocountB = occpiB[0]; 
    int vcountB = virtpiB[0];
    for(int h=1; h < nirrep_; h++) {
      occ_offB[h] = ocountB;
      ocountB += occpiB[h];
      
      vir_offB[h] = vcountB;
      vcountB += virtpiB[h];
    }

      // print
      if (print_ > 1) {
          for(int h = 0; h < nirrep_; h++) {
              fprintf(outfile," h, occ_offB[h]: %2d %2d \n", h, occ_offB[h]);
              fflush(outfile);
          }
          fprintf(outfile,"\n");
          fflush(outfile);

          for(int h = 0; h < nirrep_; h++) {
              fprintf(outfile," h, vir_offB[h]: %2d %2d \n", h, vir_offB[h]);
              fflush(outfile);
          }
          fprintf(outfile,"\n");
          fflush(outfile);
      }
 
/********************************************************************************************/
/************************** Read orbital coefficients ***************************************/
/********************************************************************************************/
        // read orbital coefficients from chkpt
        Cb_ = SharedMatrix(reference_wavefunction_->Cb());
	Cb_ref = boost::shared_ptr<Matrix>(new Matrix("Ref beta MO coefficients", nirrep_, nsopi_, nmopi_));
	
	// read orbital coefficients from external files
	if (read_mo_coeff == "TRUE"){
	  fprintf(outfile,"\n\tReading MO coefficients in pitzer order from the external file CmoB.psi...\n");  
	  fflush(outfile);
	  double **C_pitzerB = block_matrix(nso_,nmo_);
	  memset(C_pitzerB[0], 0, sizeof(double)*nso_*nmo_);

	  // read binary data
	  ifstream InFile2;
	  InFile2.open("CmoB.psi", ios::in | ios::binary);
	  InFile2.read( (char*)C_pitzerB[0], sizeof(double)*nso_*nmo_);
	  InFile2.close();
	
	  //set C_scf
	  Cb_->set(C_pitzerB);
	  free_block(C_pitzerB);
        }
        
        // Build Reference MOs
	Cb_ref->copy(Cb_);
	if(print_ > 2) Cb_->print();
	
}// end if (reference_ == "UNRESTRICTED") 

//fprintf(outfile,"\n get_moinfo is done. \n"); fflush(outfile);
}// end get_moinfo
}} // End Namespaces

