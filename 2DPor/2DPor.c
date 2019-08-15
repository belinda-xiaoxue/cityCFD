/******************************************************************/
/****CITY SCALE CFD£¨CSCFD£©Version 2: 2D basic + porous medium****/
/* Rewrite the governing equation by adding the source derived ****/
/* after coordinate trasformation & add absorbing to the top of****/
/* the domain. The buildings are treated as porous medium**********/
/*******************************************************************/

/*header*/
#include "udf.h"
/*Parameters*/
/*thermal stability*/
#define C2 -0.007 				/*-d_theta/dz*/
#define kesai -0.00014  	/* The transformation coefficient (Wang and Li,2016)*/
/* Material property */
#define C3 0.0033   	  	/* Thermal expansion coefficient*/
#define air_den 1.225  		/* Air density */
#define c_p 1006       		/* Specific heat */
/*Heat flux for urban and rural*/
#define Q_sen_city 200 		/* urban heat flux */
#define Q_sen_rural 50		/* rural heat flux */
/*damping facotr*/
#define Z_top 2490.0     	/* The height of the domain*/
#define Z_damping 1990.0 	/* Where the damping layer begins*/
/* porous */
#define Por_0 0.75       	/* poriosity 0 */
#define H_Urb 50.0       	/* height of urban*/
#define D_p 50.0         	/* characteristic size of solid paticles or cube of building */ 
#define Por_xs 100.0     	/* porous region 1 starts in x-direction*/
#define Por_xe 5000.0 		/* porous region 1 ends   in x-direction*/
#define rural_r_s 5000.0
#define rural_r_e 55000.0
/* Turbulence ralated */
#define Ce1 1.44    			/* C1-Epsilon */
#define CMU 0.09  
#define Pr_t 0.85      		/* The turbulent Prandtl-number (0.85) */
/* Constants */
#define g_acc 9.8      		/* The gravitational acceleration*/
#define pai 3.1416     		/* Pi, ratio of the circumference of a circle to its diameter */
/* Initial condition */
#define temp_air_ini 298.0/* The inital air temperautre in the simulation */

/********  UDM define *********/
DEFINE_ON_DEMAND(UDMInit)
{ 
/* Notice: Run this first to initialize your UDM*/
	Domain *domain;   /* Declare domain pointer since it is not passed a argument to DEFINE macro  */
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
			/*poriosity*/
			x_s=x[0];
			if( y <= H_Urb)
			{
				if (x_s>=Por_xs && x_s<=Por_xe) 
				{
					tmp_por = 0.75;
					C_UDMI(c,t,1) = tmp_por;
					tmp_beta=1.0;
				}
				else
				{
					tmp_por = 1.0;
					C_UDMI(c,t,1) = tmp_por;
					tmp_beta=1.0;
				}
			}else
			{
				tmp_por = 1.0;
				C_UDMI(c,t,1) = tmp_por;
				tmp_beta=1.0;
			}
			tmp_K1=pow(D_p,2)*pow(C_UDMI(c,t,1),3.0);
			tmp_K2=150*pow((1-C_UDMI(c,t,1)),2.0);	
			C_UDMI(c,t,2)=tmp_K2/tmp_K1;
			tmp_C1=D_p*pow(C_UDMI(c,t,1),3.0);
			tmp_C2=3.5*(1-C_UDMI(c,t,1))*tmp_beta/tmp_C1;
			C_UDMI(c,t,3)= tmp_C2;
		}
		end_c_loop (c,t)
	}
}
/********  profile *********/
DEFINE_PROFILE(porosity_function_Hang,t,i)
{
cell_t c;
begin_c_loop(c,t)
{
	C_PROFILE(c,t,i) = C_UDMI(c,t,1);
}
end_c_loop(c,t)
}
/****		Viscous Resistance		****/
DEFINE_PROFILE(vis_res_Hang,t,i)
{
cell_t c;
begin_c_loop(c,t)
{
	F_PROFILE(c,t,i) =  C_UDMI(c,t,2);  /*1/K*/;
}
end_c_loop(c,t)
}
/****		Initial Resistance		****/
DEFINE_PROFILE(ini_res_Hang,t,i)
{
cell_t c;
begin_c_loop(c,t)
{
	F_PROFILE(c,t,i) = C_UDMI(c,t,3); /*C2*/	
}
end_c_loop(c,t)
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
	con_por=(1-C_UDMI(c,t,1))*Q_sen_city/H_Urb;
	source = con+damping_factor*(C_T(c,t)-temp_air_ini)+con_por;
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
	tmp_K=C_UDMI(c,t,1);
	tmp_C2=C_UDMI(c,t,2);;
	Q=sqrt(u*u+v*v);
	Fk_1=u*u*dudx+v*v*dvdy+u*v*(dudy+dvdx);
	mu_t1 = rho*CMU*k/d;
	Fk=mu_t1*Fk_1/Q;
	con1=-2*Por_0*mu_l*tmp_K;
	con2= (-4/3)*pow(Por_0,2.0)*tmp_C2*rho*Q;
	con3=pow(Por_0,2.0)*tmp_C2*Fk;
	con=con1*k+con2*k+con3*k;
	source = con+con_stra;
	dS[eqn] = con1+con2+con3*2;
	
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
	tmp_K=C_UDMI(c,t,1);
	tmp_C2=C_UDMI(c,t,2);;
	Q=sqrt(u*u+v*v);
	con1=-2*Por_0*mu_l*tmp_K;
	con2= (-4/3)*pow(Por_0,2)*tmp_C2*rho*Q;
	con=con1*ed+con2*ed;
	source = con+con_stra*C_D(c,t);
	dS[eqn] = con1+con2+con_stra;

	return source;
}
/******** city sensible heat flux  *********/
DEFINE_PROFILE(WallFlux_city,thread,i)
{
	real source;
	face_t f;
	begin_f_loop(f,thread)
	{
		source=Por_0*Q_sen_city;
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
		source=Por_0*Q_sen_rural;
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


