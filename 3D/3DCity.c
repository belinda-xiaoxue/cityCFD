/********************************************************************/
/*******          CITY SCALE CFD£¨CSCFD£© Version 2: 3D (meso) ******/    
/* Rewrite the governing equation by adding the source derived after*/
/* coordinate trasformation & add absorbing to the top of the domain*/  
/********************************************************************/

/*header*/
#include "udf.h"

/*Parameters*/
/*thermal stability*/
#define C2 -0.005       	/*-d_theta/dz */
/*Also kown as negative number of the vertical gradient of the potential temperature (Wang,2015)*/
#define kesai -0.0001  		/* The transformation coefficient (Wang and Li,2016)*/
/* Material property */
#define C3 0.0033      		/* Thermal expansion coefficient*/
#define air_den 1.225  		/* Air density */
#define c_p 1006       		/* Specific heat */
/*Heat flux difference*/
#define Q_sen_city 200.0  /*urban heat flux */
#define Q_sen_rural 50.0	/*rural heat flux */
/* Turbulence ralated */
#define Ce1 1.44       		/* C1-Epsilon */
#define Pr_t 0.85      		/* The turbulent Prandtl-number (0.85) */
/* damping facotr */
#define Z_top 2480.0			/* The height of the domain*/
#define Z_damping 2000.0  /* Where the damping layer begins*/
/* Others */
#define g_acc 9.8      		/* The gravitational acceleration*/
#define pai 3.1416     		/* Pi, ratio of the circumference of a circle to its diameter */
/* Initial condition */
#define temp_air_ini 300.0/* The inital air temperautre in the simulation */

/* Define user defined memory */
DEFINE_ON_DEMAND(UDMInit)
{ 
/* Notice: Run this first to initialize your UDM*/
	Domain *domain;   /* Declare domain pointer since it is not passed a argument to DEFINE macro  */
	real  x[ND_ND];
	real z;
	cell_t c;
	Thread *t;
	domain=Get_Domain(1); /* Get the domain using Fluent utility */
   /* Loop over all cell threads in the domain */
	thread_loop_c (t,domain)
	{
		begin_c_loop (c,t)
		{
			C_CENTROID(x,c,t);
			z=x[2];
			C_UDMI(c,t,0)=1/(1+kesai*z);
		}
		end_c_loop (c,t)
	}
}

/* Modify the source terms and add the absorbing layer */
DEFINE_SOURCE(u_damping,c,t,dS,eqn) /*u*/
{
	real x[ND_ND];
	real z,con, source;
	real damping_factor;
	real z_ref_height;
	C_CENTROID(x,c,t);
	z=x[2];
	z_ref_height=(z-Z_damping)/(Z_top-Z_damping);
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
DEFINE_SOURCE(v_damping,c,t,dS,eqn) /*v*/
{
	real x[ND_ND];
	real z,con, source;
	real damping_factor;
	real z_ref_height;
	C_CENTROID(x,c,t);
	z=x[2];
	z_ref_height=(z-Z_damping)/(Z_top-Z_damping);
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
	con =damping_factor*C_V(c,t);
	source = con;
	dS[eqn] = damping_factor;
	return source;
}
DEFINE_SOURCE(w_damping,c,t,dS,eqn) /*w*/
{
	real x[ND_ND];
	real z,con,con2, source;
	real damping_factor;
	real z_ref_height;
	real J;
	C_CENTROID(x,c,t);
	z=x[2];
	z_ref_height=(z-Z_damping)/(Z_top-Z_damping);
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
	con =damping_factor*C_W(c,t);
	J=C_UDMI(c,t,0);
	con2=air_den*(J*J-1)*C3*g_acc*(C_T(c,t)-temp_air_ini);
	source = con+con2;
	dS[eqn] = damping_factor;
	return source;
}
DEFINE_SOURCE(energy_source,c,t,dS,eqn) /*energy*/
{
	real x[ND_ND];
	real z,con, source;
	real damping_factor;
	real z_ref_height;
	real J;
	C_CENTROID(x,c,t);
	z=x[2];
	z_ref_height=(z-Z_damping)/(Z_top-Z_damping);
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
	con = C2*air_den*c_p*C_W(c,t)*J;
	source = con+damping_factor*(C_T(c,t)-temp_air_ini);
	dS[eqn] = damping_factor;
	return source;
}
DEFINE_SOURCE(kinetic_source,c,t,dS,eqn) /*k*/
{
	real con, source;
	con = C3*g_acc*C_MU_T(c,t)*C2/Pr_t;
	source = con;
	dS[eqn] = 0;
	return source;
}
DEFINE_SOURCE(epsilon_source,c,t,dS,eqn) /*epsilon*/
{
	real con, source;
	real tmp_v;
	tmp_v=sqrt(C_U(c,t)*C_U(c,t)+C_V(c,t)*C_V(c,t));
	con = Ce1*tanh(fabs(C_W(c,t)/tmp_v))*C3*g_acc*C_MU_T(c,t)*C2/(Pr_t*C_K(c,t));
	source = con*C_D(c,t);
	dS[eqn] = con;
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