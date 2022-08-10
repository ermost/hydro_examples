/*
 * =====================================================================================
 *
 *       Filename:  reconstruct.hh
 *
 *    Description:  WENO-Z reconstruction
 *
 *        Version:  1.0
 *        Created:  13.06.2018 17:41:53
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#pragma once

template <typename T> inline constexpr T SQ(T x) { return x * x; };

template <bool plus_direction, bool FD, typename T = double>
struct WENO5_Branch{
public:
static inline T branch(T *WENOVAR) {

  int minus2 = 0;
  int minus1 = 1;
  int plus0 = 2;
  int plus1 = 3;
  int plus2 = 4;
  if (!plus_direction) {
    minus2 = 4;
    minus1 = 3;
    plus0 = 2;
    plus1 = 1;
    plus2 = 0;
  }

  constexpr T dFD[3]{1. / 16., 10. / 16., 5. / 16.};

  constexpr T dFV[3]{1. / 10., 3. / 5., 3. / 10.};

  T f[3];
  if (FD) {
    // Finite difference reconstruction

    f[0] = 3. / 8. * WENOVAR[minus2] - 10. / 8. * WENOVAR[minus1] +
           15. / 8. * WENOVAR[plus0];
    f[1] = -1. / 8. * WENOVAR[minus1] + 6. / 8. * WENOVAR[plus0] +
           3. / 8. * WENOVAR[plus1];
    f[2] = 3. / 8. * WENOVAR[plus0] + 6. / 8. * WENOVAR[plus1] -
           1. / 8. * WENOVAR[plus2];
  } else {
    f[0] = 1. / 3. * WENOVAR[minus2] - 7. / 6. * WENOVAR[minus1] +
           11. / 6. * WENOVAR[plus0];
    f[1] = -1. / 6. * WENOVAR[minus1] + 5. / 6. * WENOVAR[plus0] +
           1. / 3. * WENOVAR[plus1];
    f[2] = 1. / 3. * WENOVAR[plus0] + 5. / 6. * WENOVAR[plus1] -
           1. / 6. * WENOVAR[plus2];
  }

  // Smooth WENO weights: Note that these are from Del Zanna et al. 2007 (A.18)

  T beta[3];
  constexpr T beta_coeff[2]{13. / 12., 0.25};
  beta[0] = beta_coeff[0] *
                SQ(WENOVAR[minus2] + WENOVAR[plus0] - 2.0 * WENOVAR[minus1]) +
            beta_coeff[1] *
                SQ(WENOVAR[minus2] - 4. * WENOVAR[minus1] + 3 * WENOVAR[plus0]);

  beta[1] = beta_coeff[0] *
                SQ(WENOVAR[minus1] + WENOVAR[plus1] - 2.0 * WENOVAR[plus0]) +
            beta_coeff[1] * SQ(WENOVAR[minus1] - WENOVAR[plus1]);

  beta[2] = beta_coeff[0] *
                SQ(WENOVAR[plus0] + WENOVAR[plus2] - 2.0 * WENOVAR[plus1]) +
            beta_coeff[1] *
                SQ(3. * WENOVAR[plus0] - 4. * WENOVAR[plus1] + WENOVAR[plus2]);

  // Rescale epsilon
  constexpr T epsL = 1.e-42;

  // WENO-Z+: Acker et al. 2016

  const T tau_5 = std::fabs(beta[0] - beta[2]);

  const T indicator[3]{(tau_5 ) / (beta[0] + epsL),
                       (tau_5 ) / (beta[1] + epsL),
                       (tau_5 ) / (beta[2] + epsL)};

  T alpha[3]{1. + SQ(indicator[0]), 1. + SQ(indicator[1]),
             1. + SQ(indicator[2])};

  T alpha_sum = 0.;
  if (FD) {
#pragma unroll
    for (int i = 0; i < 3; ++i) {
      alpha[i] *= dFD[i];
      alpha_sum += alpha[i];
    };
  } else { // FV

#pragma unroll
    for (int i = 0; i < 3; ++i) {
      alpha[i] *= dFV[i];
      alpha_sum += alpha[i];
    };
  }

  T flux = 0.;
#pragma unroll
  for (int i = 0; i < 3; ++i) {
    flux += f[i] * alpha[i] / alpha_sum;
  };

  return flux;
};
};




template<typename T>
static inline T minmod(const T& A, const T& B){ return 0.5*(std::copysign(1,A)+std::copysign(1,B))*std::min(std::fabs(A),std::fabs(B));};
template<typename T>
static inline T median(const T& A, const T& B, const T& C){ return A + minmod(B-A,C-A);};

static constexpr double MP5_alpha = 4.;
static constexpr double MP5_eps = 1.e-10;

template <bool plus_direction, bool FD, typename T = double, bool fifth=true>
struct MP5_Branch{
public:
static inline T branch(T *U) {

  int minus2 = 0;
  int minus1 = 1;
  int plus0 = 2;
  int plus1 = 3;
  int plus2 = 4;
  if (!plus_direction) {
    minus2 = 4;
    minus1 = 3;
    plus0 = 2;
    plus1 = 1;
    plus2 = 0;
  }


     T Ur;
      
     if(fifth){
	// Originial Huynh & Suresh
	if(!FD){
	  Ur = (2.*U[minus2]- 13.*U[minus1]+ 47.*U[plus0]+ 27. * U[plus1] - 3.*U[plus2])/60.;
	}else{
	// Del Zanna et al. ECHO code for DER 
	  Ur = (3./128.*U[minus2]- 20./128.*U[minus1]+ 90./128.*U[plus0]+ 60./128. * U[plus1] - 5./128.*U[plus2]);
	};
     }else{
	    if(!FD){
      		  Ur = (7.0/12.0)*(U[plus1] + U[plus0] ) - (1.0/12.0)*(U[minus1]  + U[plus2]);
	    }else{
      		  Ur = (9.0/16.0)*(U[plus1] + U[plus0] ) - (1.0/16.0)*(U[minus1]  + U[plus2]);
	    };
     }


      //Eq. 2.12 in Suresh & Huynh
      T Ur_MP= U[plus0] + minmod(U[plus1]-U[plus0], MP5_alpha*(U[plus0]-U[minus1]));

      T Unorm=0.0;
      for(int i=minus2; i<=plus2 ; ++i) Unorm+=U[i]*U[i];


      //Eq. 2.30 in Suresh & Huynh
      if((Ur-U[plus0])*(Ur-Ur_MP)> MP5_eps*Unorm){
      		T d[3];
		for(int i=0; i<3;++i) d[i]= U[minus2 + i] + U[plus0 +i]- 2.*U[minus1 + i];
//		d[1]= U[minus1] + U[plus1]- 2.*U[plus0];
//		d[2]= U[plus0] + U[plus2]- 2.*U[plus1];

		const T dM4_r= minmod(minmod(4.*d[1]-d[2], 4.*d[2]-d[1]), minmod(d[1],d[2]));
		const T dM4_l= minmod(minmod(4.*d[1]-d[0], 4.*d[0]-d[1]), minmod(d[1],d[0]));

		const T deltaU= U[plus0]-U[minus1];		
		
		//Eq. 2.8
		const T U_UL= U[plus0] + MP5_alpha*deltaU;
		//Eq. 2.16
		const T U_AV= 0.5*(U[plus0] + U[plus1]);
	
		//Eq. 2.28
		const T U_MD= U_AV - 0.5*dM4_r;
		//Eq 2.29
		const T U_LC= U[plus0] + 0.5*deltaU + 4./3.*dM4_l;

		//Eq. 2.24
		const T U_min=std::max(std::min(U[plus0], std::min(U[plus1], U_MD)), std::min(U[plus0], std::max(U_UL, U_LC)));
		const T U_max=std::min(std::max(U[plus0], std::max(U[plus1], U_MD)), std::max(U[plus0], std::max(U_UL, U_LC)));

		Ur =  median(U_min, Ur, U_max);	
      }; 

      return Ur;
};
};



template <bool FD , typename T , template<bool,bool,typename> class Rec>
struct RecInterface{
  public:
    static inline void reconstruct( T* U, T& Ur, T& Ul){
	
      Ur = Rec<true, FD,T>::branch(U);
      Ul = Rec<false, FD,T>::branch(U);

    };
};


static constexpr double xppmC = 1.25;

template <bool FD = false, typename T = double>
struct XPPM_Branch{
public:
static inline void reconstruct(T *U, T& Ur, T& Ul) {

  int minus2 = 0;
  int minus1 = 1;
  int plus0 = 2;
  int plus1 = 3;
  int plus2 = 4;

     if(!FD){ 
       Ur = (7.0/12.0)*(U[plus1] + U[plus0] ) - (1.0/12.0)*(U[minus1]  + U[plus2]);
       Ul = (7.0/12.0)*(U[minus1] + U[plus0] ) - (1.0/12.0)*(U[plus1]  + U[minus2]);
     }else{
       Ur = (9.0/16.0)*(U[plus1] + U[plus0] ) - (1.0/16.0)*(U[minus1]  + U[plus2]);
       Ul = (9.0/16.0)*(U[minus1] + U[plus0] ) - (1.0/16.0)*(U[plus1]  + U[minus2]);
     }


     auto monotonize = [&]  (const int i, auto &Uf){

       //Check monotonicity!
       if( (U[i] - Uf) * (U[i+1] - Uf) > 0) {
	    T D2U = 3.*(U[i] - 2.*Uf + U[i+1]);
	    T D2UL = (U[i-1] - 2.*U[i] + U[i+1]);
	    T D2UR = (U[i] - 2.*U[i+1] + U[i+2]);

            T D2Ulim = 0.;
	    if((std::signbit(D2U) ==  std::signbit(D2UL)) && (std::signbit(D2UL) ==  std::signbit(D2UR))){
	      D2Ulim = std::copysign(1.,D2U)* std::min(std::fabs(D2U) , xppmC * std::min(std::fabs(D2UL),std::fabs(D2UR)));	    
	    }
	    
	    Uf = 0.5*(U[i] + U[i+1]) - 1./6.*D2Ulim;
       }
     };

     monotonize(plus0, Ur);
     monotonize(minus1, Ul);



     //Check for extrema!
     if( ((Ur - U[plus0])*(U[plus0] - Ul) <= 0.)
         || (U[minus1] - U[plus0])*(U[plus0] - U[plus1]) <=0 ){

       T D2U = 6.*(Ul - 2.*U[plus0] + Ur);
       T D2UC = (U[minus1] - 2.*U[plus0] + U[plus1]);
       T D2UL= (U[minus2] - 2.*U[minus1] + U[plus0]);
       T D2UR= (U[plus0] - 2.*U[plus1] + U[plus2]);

       T D2Ulim =0;
       if((std::signbit(D2U) ==  std::signbit(D2UL)) && 
	   (std::signbit(D2UL) ==  std::signbit(D2UR) &&
	   (std::signbit(D2UC) ==  std::signbit(D2UR)
	     ))){
	      D2Ulim = std::copysign(1.,D2U)* std::min(std::fabs(D2U) , xppmC * std::min(std::fabs(D2UL),std::min(std::fabs(D2UR),std::fabs(D2UC))));	    
 	};
	if(D2U !=0){
		auto frac  =D2Ulim/D2U;
//		frac = std::copysign(1.,frac)*std::max(1.e-12, std::fabs(frac));
		Ur = U[plus0] + (Ur - U[plus0])*frac;
		Ul = U[plus0] + (Ul - U[plus0])*frac;
	}
     } else{

       //Limit variations
       //
      
       auto deltaR = Ur - U[plus0];
       auto deltaL = Ul - U[plus0];
     
       if(std::fabs(deltaR) >= 3.* std::fabs(deltaL)){
	   Ur = U[plus0] - 3.* deltaL;
       };

       if(std::fabs(deltaL) >= 3.* std::fabs(deltaR)){
	   Ul = U[plus0] - 3.* deltaR;
       };
     };
   };
};



template <bool FD, typename T , template<bool,typename> class Rec> class FifthOrder_Reconstruct {

public:
  template <int ndim, int dir, typename Tstorage>
  static inline void reconstruct(Tstorage &U, Tstorage &scratchR,
                                 Tstorage &scratchL) {

    if (ndim == 1) {
//      assert(dir==0);

      for (int nv = 0; nv < Tstorage::nsystem; ++nv){
#pragma omp simd
        for (int ijk = 2; ijk < U.grid.extent[0] - 2; ++ijk) {
          const int i = ijk + nv * U.grid.extent[0];
          T Ulocal[5]{U.U[i - 2], U.U[i - 1], U.U[i - 0], U.U[i + 1],
                      U.U[i + 2]};

          // The shift is because we want L and R to be at the same interface
          // for the same index
//          scratchR.U[i] = branch<false, FD, T>::branch(Ulocal);
//          scratchL.U[i+1] = branch<true, FD, T>::branch(Ulocal);

	    Rec<FD,T>::reconstruct(Ulocal, scratchL.U[i+1], scratchR.U[i]);	

	};
      };
    };

    if (ndim == 2) {

      if (dir == 0) {

        for (int nv = 0; nv < Tstorage::nsystem; ++nv)
#pragma omp parallel for
          for (int j = 0; j < U.grid.extent[1]; ++j)
            for (int i = 2; i < U.grid.extent[0] - 2; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);

              T Ulocal[5]{U.U[ij - 2], U.U[ij - 1], U.U[ij - 0], U.U[ij + 1],
                          U.U[ij + 2]};

              // The shift is because we want L and R to be at the same
              // interface for the same index
//              scratchR.U[ij] = branch<false, FD, T>::branch(Ulocal);
//              scratchL.U[ij+1] = branch<true, FD, T>::branch(Ulocal);
	    Rec<FD,T>::reconstruct(Ulocal, scratchL.U[ij+1], scratchR.U[ij]);	
            };

      }; // dir ==0

      if (dir == 1) {

        for (int nv = 0; nv < Tstorage::nsystem; ++nv)
#pragma omp parallel for
          for (int j = 2; j < U.grid.extent[1] - 2; ++j)
            for (int i = 0; i < U.grid.extent[0]; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int jO = U.grid.extent[0];

              T Ulocal[5]{U.U[ij - 2 * jO], U.U[ij - 1 * jO], U.U[ij - 0],
                          U.U[ij + 1 * jO], U.U[ij + 2 * jO]};

              // The shift is because we want L and R to be at the same
              // interface for the same index
//              scratchR.U[ij ] = branch<false, FD, T>::branch(Ulocal);
//              scratchL.U[ij+1*jO] = branch<true, FD, T>::branch(Ulocal);
	    Rec<FD,T>::reconstruct(Ulocal, scratchL.U[ij+1*jO], scratchR.U[ij]);	
            };

      }; // dir ==1

    }; // ndim==2

    if (ndim == 3) {
      assert(!"Not implemented yet");
    }

    //Some fixing 
      scratchR.is_primitive = U.is_primitive;
      scratchL.is_primitive = U.is_primitive;


  };
};


template <bool FD = false, typename T = double> 
using IFWENO5 = RecInterface<FD,T, WENO5_Branch>;

template <bool FD = false, typename T = double> 
using WenoZ_Reconstruct =  FifthOrder_Reconstruct<FD,T,IFWENO5>;

template <bool dir, bool FD , typename T = double> 
using MP5_Branch5 = MP5_Branch<dir,FD,T,true>;

template <bool dir, bool FD , typename T = double> 
using MP5_Branch4 = MP5_Branch<dir,FD,T,false>;

template <bool FD = false, typename T = double> 
using IFMP5 = RecInterface<FD,double,MP5_Branch5>;

template <bool FD = false, typename T = double> 
using IFMP54 = RecInterface<FD,double,MP5_Branch4>;

template <bool FD = false, typename T = double> 
using MP5_Reconstruct =  FifthOrder_Reconstruct<FD,T, IFMP5>;

template <bool FD = false, typename T = double> 
using MP54_Reconstruct =  FifthOrder_Reconstruct<FD,T, IFMP54>;

template <bool FD = false, typename T = double> 
using XPPM_Reconstruct =  FifthOrder_Reconstruct<FD,T, XPPM_Branch>;


template <typename Rec,int Boffset=0, bool FD = false, typename T = double> class ReconstructMHD_Special {

public:
  template <int ndim, int dir, typename Tstorage>
  static inline void reconstruct(Tstorage &U, Tstorage &scratchR,
                                 Tstorage &scratchL) {

    Rec::template reconstruct<ndim,dir>(U,scratchR,scratchL);

    auto const rec_mhd = [&] (int const nv){

    if (ndim == 1) {
//      assert(dir==0);

#pragma omp simd
        for (int ijk = 2; ijk < U.grid.extent[0] - 2; ++ijk) {
          const int i = ijk + nv * U.grid.extent[0];
          T Ulocal[5]{U.U[i - 1], U.U[i - 0], U.U[i + 1],
                      U.U[i + 2]};

          // The shift is because we want L and R to be at the same interface
          // for the same index
	  if(!FD){
          	scratchR.U[i] = -1./12.*U[i-2] + 7./12.*U[i-1] + 7./12.*U[i] - 1./12. * U[i+1];
	  }else{
          	scratchR.U[i] = -1./16.*U[i-2] + 9./16.*U[i-1] + 9./16.*U[i] - 1./16. * U[i+1];
	  };
          scratchL.U[i] = scratchR.U[i];

	};
    };

    if (ndim == 2) {

      if (dir == 0) {

#pragma omp parallel for
          for (int j = 0; j < U.grid.extent[1]; ++j)
            for (int i = 2; i < U.grid.extent[0] - 2; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * (nv));

              T Ulocal[4]{U.U[ij - 2], U.U[ij - 1], U.U[ij - 0], U.U[ij + 1]};

	    // The shift is because we want L and R to be at the same interface
	    // for the same index
	    if(!FD){
		  scratchR.U[i] = -1./12.*Ulocal[0] + 7./12.*Ulocal[1] + 7./12.*Ulocal[2] - 1./12. * Ulocal[3];
	    }else{
		  scratchR.U[i] = -1./16.*Ulocal[0] + 9./16.*Ulocal[1] + 9./16.*Ulocal[2] - 1./16. * Ulocal[3];
	    };
	    scratchL.U[i] = scratchR.U[i];
            };

      }; // dir ==0

      if (dir == 1) {

#pragma omp parallel for
          for (int j = 2; j < U.grid.extent[1] - 2; ++j)
            for (int i = 0; i < U.grid.extent[0]; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * (nv));
              const int jO = U.grid.extent[0];

              T Ulocal[4]{U.U[ij - 2 * jO], U.U[ij - 1 * jO], U.U[ij - 0], U.U[ij + 1 * jO]};

	      // The shift is because we want L and R to be at the same interface
	      // for the same index
	      if(!FD){
		    scratchR.U[i] = -1./12.*Ulocal[0] + 7./12.*Ulocal[1] + 7./12.*Ulocal[2] - 1./12. * Ulocal[3];
	      }else{
		    scratchR.U[i] = -1./16.*Ulocal[0] + 9./16.*Ulocal[1] + 9./16.*Ulocal[2] - 1./16. * Ulocal[3];
	      };
	      scratchL.U[i] = scratchR.U[i];
	      };

      }; // dir ==1

    }; // ndim==2

    if (ndim == 3) {
      assert(!"Not implemented yet");
    }
    };

    rec_mhd(Boffset+dir);
    rec_mhd(Boffset+3); //Div cleaning

    //Some fixing 
      scratchR.is_primitive = U.is_primitive;
      scratchL.is_primitive = U.is_primitive;


  };
};

template <typename T>
struct TVD_MinMod{
  public:
    static inline void reconstruct( T* U, T& Ur, T& Ul){

      auto deltap = U[2] - U[1];
      auto deltam = -(U[0] - U[1]);
	
      Ur = U[1] + 0.5*minmod(deltam,deltap);
      Ul = U[1] + 0.5*minmod(-deltap,-deltam);

    };
};

template <bool FD , typename T , template<typename> class Rec> class Reconstruct_TVD {

public:
  template <int ndim, int dir, typename Tstorage>
  static inline void reconstruct(Tstorage &U, Tstorage &scratchR,
                                 Tstorage &scratchL) {

    if (ndim == 1) {
//      assert(dir==0);

      for (int nv = 0; nv < Tstorage::nsystem; ++nv){
#pragma omp simd
        for (int ijk = 1; ijk < U.grid.extent[0] - 1; ++ijk) {
          const int i = ijk + nv * U.grid.extent[0];
          T Ulocal[3]{U.U[i - 1], U.U[i - 0], U.U[i + 1]};

	    Rec<T>::reconstruct(Ulocal, scratchL.U[i+1], scratchR.U[i]);	

	};
      };
    };

    if (ndim == 2) {

      if (dir == 0) {

        for (int nv = 0; nv < Tstorage::nsystem; ++nv)
#pragma omp parallel for
          for (int j = 0; j < U.grid.extent[1]; ++j)
            for (int i = 1; i < U.grid.extent[0] - 1; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);

              T Ulocal[3]{U.U[ij - 1], U.U[ij - 0], U.U[ij + 1]};

              // The shift is because we want L and R to be at the same
              // interface for the same index
	    Rec<T>::reconstruct(Ulocal, scratchL.U[ij+1], scratchR.U[ij]);	
            };

      }; // dir ==0

      if (dir == 1) {

        for (int nv = 0; nv < Tstorage::nsystem; ++nv)
#pragma omp parallel for
          for (int j = 1; j < U.grid.extent[1] - 1; ++j)
            for (int i = 0; i < U.grid.extent[0]; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int jO = U.grid.extent[0];

              T Ulocal[3]{U.U[ij - 1 * jO], U.U[ij - 0], U.U[ij + 1 * jO]};

              // The shift is because we want L and R to be at the same
              // interface for the same index
	    Rec<T>::reconstruct(Ulocal, scratchL.U[ij+1*jO], scratchR.U[ij]);	
            };

      }; // dir ==1

    }; // ndim==2

    if (ndim == 3) {
      assert(!"Not implemented yet");
    }

    //Some fixing 
      scratchR.is_primitive = U.is_primitive;
      scratchL.is_primitive = U.is_primitive;


  };
};

template <bool FD = false, typename T = double> 
using MinMod_Reconstruct =  Reconstruct_TVD<FD,T, TVD_MinMod>;



