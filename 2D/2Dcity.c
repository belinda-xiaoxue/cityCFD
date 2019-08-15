/********************************************************************/
/*******          CITY SCALE CFD£¨CSCFD£© Version 1: 2D basic  ******/    
/* Rewrite the governing equation by adding the source derived after*/
/* coordinate trasformation & add absorbing to the top of the domain*/  
/********************************************************************/

/*header*/
#include "udf.h"

/*Parameters*/
/*thermal stability*/
#define C2 -0.007 				/*-d_theta/dz*/
/*Also kown as negative number of the vertical gradient of the potential temperature (Wang,2015)*/
#define kesai -0.00014    /* The transformation coefficient (Wang and Li,2016)*/
/* Material property */
#define C3 0.0033  /* Thermal expansion coefficient*/
#define air_den 1.225  		/* Air density */
#define c_p 1006       		/* Specific heat */
/*Heat flux for urban and rural*/
#define Q_sen_city 200    /* urban heat flux */
#define Q_sen_rural 50 		/* rural heat flux */
/*damping facotr*/
#define Z_top 2490.0 			/* The height of the domain*/
#define Z_damping 1990.0	/* Where the damping layer begins*/
/* Turbulence ralated */
#define Ce1 1.44    			/* C1-Epsilon */
#define CMU 0.09  
#define Pr_t 0.85      		/* The turbulent Prandtl-number (0.85) (Kristof et al.,2009) */
/* Constants */
#define g_acc 9.8      		/* The gravitational acceleration*/
#define pai 3.1416     		/* Pi, ratio of the circumference of a circle to its diameter */
/* Initial condition */
#define temp_air_ini 298.0/* The inital air temperautre in the simulation */
/********  UDM define *********/
DEFINE_ON_DEMAND(UDMInit)
{ 
/* run this first to initialize your UDM*/
	Domain *domain;   /* declare domain pointer since it is not passed a argument to DEFINE macro  */
	real  x[ND_ND];
	real y,x_s,tmp_por;
	real tmp_K1,tmp_K2,tmp_C1,tmp_C2,tmp_beta;
	cell_t c;
	Thread *t;
	domain=Get_Domain(1); /* Get the domain using Fluent utility */
   /* Loop over all cell threads in the domain */
	thread_loop_c (t,domain)
	{
		begin_c_loop (c,t)
		{
			C_CENTROID(x,c,t);
			/*vertical cooridnate parameter*/
			y=x[1];
			C_UDMI(c,t,0)=1/(1+kesai*y);
		}
		end_c_loop (c,t)
	}
}

/********  source term *********/
DEFINE_SOURCE(u_damping,c,t,dS,eqn)
{
	real x[ND_ND];
	real y,con, source;
	real damping_factor;
	real z_ref_height;
	C_CENTROID(x,c,t);
	y=x[1];
	z_ref_height=(y-Z_damping)/(Z_top-Z_damping);
	if (z_ref_height<=0)
	{
		damping_factor=0;
	}
	else if (z_ref_height<=0.5 && z_ref_height>0)
	{
		damping_factor=-1*(1-cos(z_ref_height*pai))/2.0;
	}
	else 
	{
		damping_factor=-1*(1+(z_ref_height-0.5)*pai)/2.0;
	}
	con =damping_factor*C_U(c,t);
	source = con;
	dS[eqn] = damping_factor;
	return source;
}
DEFINE_SOURCE(v_damping,c,t,dS,eqn)
{
	real x[ND_ND];
	real y,con,con2,source;
	real damping_factor;
	real z_ref_height;
	real J;
	C_CENTROID(x,c,t);
	y=x[1];
	z_ref_height=(y-Z_damping)/(Z_top-Z_damping);
	if (z_ref_height<=0)
	{
		damping_factor=0;
	}
	else if (z_ref_height<=0.5 &&z_ref_height>0.5)
	{
		damping_factor=-1*(1-cos(z_ref_height*pai))/2.0;
	}
	else 
	{
		damping_factor=-1*(1+(z_ref_height-0.5)*pai)/2.0;
	}
	con =damping_factor*C_V(c,t);
	J=C_UDMI(c,t,0);
	con2=air_den*(J*J-1)*C3*g_acc*(C_T(c,t)-temp_air_ini);
	source = con+con2;
	dS[eqn] = damping_factor;
	return source;
}
DEFINE_SOURCE(energy_source,c,t,dS,eqn)
{
	real x[ND_ND];
	real y,con, source;
	real damping_factor;
	real z_ref_height;
	real con_por;
	real J;
	C_CENTROID(x,c,t);
	y=x[1];
	z_ref_height=(y-Z_damping)/(Z_top-Z_damping);
	if (z_ref_height<=0)
	{
		damping_factor=0;
	}
	else if (z_ref_height<=0.5 && z_ref_height>0)
	{
		damping_factor=-1*(1-cos(z_ref_height*pai))/2.0;
	}
	else 
	{
		damping_factor=-1*(1+(z_ref_height-0.5)*pai)/2.0;
	}
	J=C_UDMI(c,t,0);
	con = C2*air_den*c_p*C_V(c,t)*J;
	source = con+damping_factor*(C_T(c,t)-temp_air_ini);
	dS[eqn] = damping_factor;
	return source;
}
/****		TKE   source				****/
DEFINE_SOURCE(kinetic_source,c,t,dS,eqn)
{
	real source;
	real tmp_K, tmp_C2;
	real con,con1,con2,con3,con_stra;
	real Fk,Q,mu_t1,Fk_1;
	real rho = C_R(c,t);
	real mu_l=C_MU_L(c,t);
	real k = C_K(c,t);
	real d = C_D(c,t);
	real u = C_U(c,t);
	real v = C_V(c,t);
	real dudx =C_DUDX(c,t);
	real dudy =C_DUDY(c,t);
	real dvdx =C_DVDX(c,t);
    real dvdy =C_DVDY(c,t);

	con_stra = C3*g_acc*C_MU_T(c,t)*C2/Pr_t;
	source = con_stra;
	dS[eqn] = 0;
	
	return source;
}
/****		TDR   source				****/
DEFINE_SOURCE(disspation_source,c,t,dS,eqn)
{
	real source;
	real tmp_K, tmp_C2;
	real con,con1,con2,con_stra;
	real Q;
	real rho = C_R(c,t);
	real mu_l=C_MU_L(c,t);
	real ed=C_D(c,t);
	real u = C_U(c,t);
	real v= C_V(c,t);

	con_stra = Ce1*tanh(C_V(c,t)/C_U(c,t))*C3*g_acc*C_MU_T(c,t)*C2/(Pr_t*C_K(c,t));
	source = con_stra*C_D(c,t);
	dS[eqn] = con_stra;

	return source;
}
/******** city sensible heat flux  *********/
DEFINE_PROFILE(WallFlux_city,thread,i)
{
	real source;
	face_t f;
	begin_f_loop(f,thread)
	{
		source=Q_sen_city;
		F_PROFILE(f,thread,i) = source;
	}
	end_f_loop(f,thread)
}
/********  rural sensible heat flux *********/
DEFINE_PROFILE(WallFlux_rural,thread,i)
{
	real source;
	face_t f;
	begin_f_loop(f,thread)
	{
		source=Q_sen_rural;
		F_PROFILE(f,thread,i) = source;
	}
	end_f_loop(f,thread)
}
DEFINE_PROFILE(WallFlux_build,thread,i)
{	
	real source;
	face_t f;
	
	begin_f_loop(f,thread)
	{
		source=Q_sen_city;
		F_PROFILE(f,thread,i) = source;
	}
	end_f_loop(f,thread)

}


