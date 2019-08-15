/******************************************************************/
/****   CITY SCALE CFD£¨CSCFD£©Version 3: 2D UHIC evolution    ****/
/* The interested building region is fully meshed and the other   */
/* city areas are treated as porous medium. Consider a UHIC       */
/* evolution for day time   (ref. Wang and Li, 2017)              */
/******************************************************************/

/*header*/
#include "udf.h"
#include "stdio.h"
#include "math.h"
#include "mem.h"
/*parameters*/
/*thermal stability*/
#define C2 -0.007						/*-d_theta/dz*/
#define kesai -0.00014			/* The transformation coefficient (Wang and Li,2016)*/
/*damping facotr*/
#define Z_top 2540.0
#define Z_damping 2040.0 
/***************************/
#define Q_anth_max 100     	/* anthropogenic heat */
/*********city area*****************/
/* porous */
#define Por_0 0.75     			/* poriosity 0 */
#define H_Urb 50.0     			/* height of urban*/
#define D_p 50.0    				/* characteristic size of solid paticles or cube of building */ 
#define Por_xs 0.0     			/* porous region 1 starts in x-direction*/
#define Por_xe 2400.0  			/* porous region 1 ends   in x-direction*/
#define Por_xs2 2600.0 			/* porous region 2 starts in x-direction*/
#define Por_xe2 5000.0 			/* porous region 2 ends   in x-direction*/
#define rural_r_s 5000.0
#define rural_r_e 55000.0
/*lines*/
#define line_2m 2.0  				/*line for y=2m*/
#define line_10m 10.0 			/*line for y=10m*/
#define eps_2m 0.7					/* the size of the mesh */
#define eps_10m 2.0					/* the size of the mesh */
/*solar*/
#define start_hour 6.5
#define updata_time_tl 600.0
#define iter_factor 1.
#define Alpha_atm 0.1		 
#define C_A 0.1			  
/*long wave*/
#define Sky_RH 0.25 				/*RH*/
/*sensible*/
#define Temp_air_init 298.0 /*initial air temperature*/
#define Temp_land_ini 293.0	/*iinitial surface temperature*/
/*surface heat flux*/
/*solar*/
#define Alpha_land_c 0.18	 	/*albedo of the surface in the city*/ 
#define Sky_VF_c 0.3				/*sky view factor in the city*/
#define Sky_emi_c 0.88 			/*emissivity in the city*/
/*sensible*/
/*#define U_ref_c 2.0				    Refernece height at 10m*/
/*soil*/
#define Deep_land_c 0.08		/**/
#define CP_land_c 500				/**/
/*********rural area*****************/
/*surface heat flux*/
/*solar*/
#define Alpha_land_r 0.4		/*albedo of the surface in the city*/ 
#define Sky_VF_r 1.0  			/*sky view factor in the city*/
#define Sky_emi_r 0.98 			/*emissivity in the city*/
/*sensible*/
/*#define U_ref_r 4.0				     Refernece height at 10m*/
/*soil*/
#define Deep_land_r 0.08		/**/
#define CP_land_r 350				/**/
/****************************/
/*constants*/
/*thermal stability*/
#define C3 0.0033
#define Ce1 1.44
#define CMU 0.09  
#define Pr_t 0.85      		/* The turbulent Prandtl-number (0.85) */
/*surface heat flux*/
/*solar*/
#define Solar_day 67
#define Soalr_lati 22.267
#define Solar_const 1368 		/**/
/*long wave*/
#define Blotz 0.0000000567
#define Dew_a 17.27 				/**/
#define Dew_b 237.7 				/**/
/*sensibel*/
#define Den_air 1.225				/**/
#define CP_air 1006.4 			/**/
#define Von_karman 0.41 		/**/
/*others*/
#define secs 3600.0
#define pi 3.14159
#define g_acc 9.8      			/* The gravitational acceleration*/

/**************************************/
/*define global variables for 2m air temperature*/
	real Q_value_sen_city=0.;
	real Q_value_sen_rural=0.;
	real Q_value_sen_build=0.;
	real ave_temp_2m_c=Temp_air_init;
	real ave_temp_2m_r=Temp_air_init;
	real surf_temp_r=Temp_land_ini;
	real surf_temp_c=Temp_land_ini;
	real surf_temp_p=Temp_land_ini;
	real U_ref_c=0.0;
	real U_ref_r=0.0;
/*******************************/
/****************Before calculation*******************************/
/************************in each calculation****************/
/*------------------------------------------------------------*/
/* -------------FUNCTION: average among the 2m line-------- */
/* ---Returns a average temperature of the whole line-------- */
/*------------------------------------------------------------*/
	double average_temp_zone(real line_meter,real eps_value, real x_start,real x_end)
	{
		int surface_thread_id=2;
		real total_volume=0.0;
		real total_temp=0.0;
		real aver_temp=0.0;
	#if !RP_HOST 
		cell_t c;
		Thread * tc;
		real xc[ND_ND];
		Domain *domain;
		domain=Get_Domain(1);
	#endif /* !RP_HOST */
		host_to_node_int_1(surface_thread_id);
   #if !RP_HOST 
		tc = Lookup_Thread(domain,surface_thread_id);
		begin_c_loop_int (c, tc)
		{
			C_CENTROID(xc,c,tc);
			if (xc[0]>=x_start && xc[0]<=x_end)
				{
				if (fabs(xc[1]-line_meter)<eps_value) 
				{
				total_volume += C_VOLUME(c, tc);
				total_temp += C_T(c, tc)*C_VOLUME(c, tc);
				}
			}
		}
		end_c_loop_int (c, tc)
	# if RP_NODE /* Perform node synchronized actions here
	Does nothing in Serial */
	total_volume = PRF_GRSUM1(total_volume);
	total_temp = PRF_GRSUM1(total_temp);
	aver_temp=total_temp/total_volume;
	# endif /* RP_NODE */
	#endif /* !RP_HOST */
	node_to_host_real_1(aver_temp); /* Does nothing in SERIAL */
	
	#if !RP_NODE  /* SERIAL or HOST */
	/*Message("\nAverage air temperature  is %f (K)\n",aver_temp);*/
	#endif /* !RP_NODE */
	return aver_temp;
	}
/*------------------------------------------------------------*/
/* -------------FUNCTION: average among the line-------- */
/* --Returns a average horizontal velocity of the whole line--- */
/*------------------------------------------------------------*/
double average_abs_vel_zone(real line_meter,real eps_value, real x_start,real x_end)
	{
		int surface_thread_id=2;
		real total_volume=0.0;
		real total_temp=0.0;
		real aver_temp=0.0;
	#if !RP_HOST 
		cell_t c;
		Thread * tc;
		real xc[ND_ND];
		Domain *domain;
		domain=Get_Domain(1);
	#endif /* !RP_HOST */
		host_to_node_int_1(surface_thread_id);
   #if !RP_HOST 
		tc = Lookup_Thread(domain,surface_thread_id);
		begin_c_loop_int (c, tc)
		{
			C_CENTROID(xc,c,tc);
			if (xc[0]>=x_start && xc[0]<=x_end)
				{
				if (fabs(xc[1]-line_meter)<eps_value) 
				{
				total_volume += C_VOLUME(c, tc);
				total_temp += C_U(c, tc)*C_VOLUME(c, tc);
				}
			}
		}
		end_c_loop_int (c, tc)
	# if RP_NODE /* Perform node synchronized actions here
	Does nothing in Serial */
	total_volume = PRF_GRSUM1(total_volume);
	total_temp = PRF_GRSUM1(total_temp);
	aver_temp=fabs(total_temp/total_volume);
	# endif /* RP_NODE */
	#endif /* !RP_HOST */
	node_to_host_real_1(aver_temp); /* Does nothing in SERIAL */
	
	#if !RP_NODE  /* SERIAL or HOST */
	/*Message("\nAverage velocity  is %f (K)\n",aver_temp);*/
	#endif /* !RP_NODE */
	return aver_temp;
	}
/*------------------------------------------------------------*/
/* -------FUNCTION: calculate the anthropogetic heat------ */
/* -------------Returns the tempearture of the surface-------- */
/*------------------------------------------------------------*/
double anthro_heat(real calc_time)
{
	real Q_anth;
	#if !RP_HOST 
		if (calc_time<=1.*secs)	  {Q_anth=0.19*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=2.*secs&&calc_time>1.*secs)	  {Q_anth=0.16*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=3.*secs&&calc_time>2.*secs)	  {Q_anth=0.08*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=4.*secs&&calc_time>3.*secs)	  {Q_anth=0.09*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=5.*secs&&calc_time>4.*secs)	  {Q_anth=0.09*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=6.*secs&&calc_time>5.*secs)	  {Q_anth=0.21*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=7.*secs&&calc_time>6.*secs)	  {Q_anth=0.45*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=8.*secs&&calc_time>7.*secs)	  {Q_anth=0.86*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=9.*secs&&calc_time>8.*secs)	  {Q_anth=0.99*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=10.*secs&&calc_time>9.*secs)	{Q_anth=0.88*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=11.*secs&&calc_time>10.*secs)	{Q_anth=0.97*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=12.*secs&&calc_time>11.*secs)	{Q_anth=0.87*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=13.*secs&&calc_time>12.*secs)	{Q_anth=0.85*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=14.*secs&&calc_time>13.*secs)	{Q_anth=0.85*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=15.*secs&&calc_time>14.*secs)	{Q_anth=0.82*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=16.*secs&&calc_time>15.*secs)	{Q_anth=0.85*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=17.*secs&&calc_time>16.*secs)	{Q_anth=0.95*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=18.*secs&&calc_time>17.*secs)	{Q_anth=0.95*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=19.*secs&&calc_time>18.*secs)	{Q_anth=0.88*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=20.*secs&&calc_time>19.*secs)	{Q_anth=0.98*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=21.*secs&&calc_time>20.*secs)	{Q_anth=0.96*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=22.*secs&&calc_time>21.*secs)	{Q_anth=0.79*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=23.*secs&&calc_time>22.*secs)	{Q_anth=0.58*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=24.*secs&&calc_time>23.*secs)	{Q_anth=0.33*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=25.*secs&&calc_time>24.*secs)	  {Q_anth=0.19*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=26.*secs&&calc_time>25.*secs)	  {Q_anth=0.16*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=27.*secs&&calc_time>26.*secs)	  {Q_anth=0.08*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=28.*secs&&calc_time>27.*secs)	  {Q_anth=0.09*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=29.*secs&&calc_time>28.*secs)	  {Q_anth=0.09*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=30.*secs&&calc_time>29.*secs)	  {Q_anth=0.21*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=31.*secs&&calc_time>30.*secs)	  {Q_anth=0.45*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=32.*secs&&calc_time>31.*secs)	  {Q_anth=0.86*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
		else if (calc_time<=33.*secs&&calc_time>32.*secs)	  {Q_anth=0.99*Q_anth_max;}/*Message("\nThe anthropogetic heat at time %g is %g\n",calc_time,Q_anth);*/
	#endif /* !RP_HOST */
return Q_anth;

}
/*------------------------------------------------------------*/
/* -------FUNCTION: calculate the surface temerature------ */
/* -------------Returns the tempearture of the surface-------- */
/*------------------------------------------------------------*/
double surface_temp(int indicator_city, real ave_temp_2m, real U_ref)
{
	real beta_tmp,Dew_temp,e_at,B_e,B_lw;
	real SW;
	real Co_lw,Co_up,Co_cons_1,Co_cons_2,Co_sh;
	real tmp_exp;
	real time, time_start;
	real tmp_r,tmp_n,t_r,temp_ref_1,temp_ref_2;
	real tmp_exp_1,tmp_exp_2;
	real x,y;	
	real Co_solar_sw1;
	real omiga,angle_lati,angle_dec;
	real time_rise,time_set,time_rise2;
	real Co_solar_sw_sin,Co_solar_sw_cos,Co_soil;
	real surf_temp=0.0;
#if !RP_HOST  /* SERIAL or NODE */
/*calculate the surface temperature*/
	omiga=pi/(12*secs);
	angle_lati=Soalr_lati/180*pi;
	angle_dec=23.45*sin(2*pi*(Solar_day-81)/365)/180*pi;
	time_rise=(acos(tan(angle_lati)*tan(angle_dec))/omiga);
	time_set=(24-time_rise/3600)*3600;
	time_rise2=time_rise+24*3600;

	time_start=start_hour*secs;
	time = CURRENT_TIME*iter_factor+time_start;	
	beta_tmp=Dew_a*(ave_temp_2m-273.15)/(Dew_b+ave_temp_2m-273.15)+log(Sky_RH);
	Dew_temp=Dew_b*beta_tmp/(Dew_a-beta_tmp);
	e_at=0.004*Dew_temp+0.8;
	B_e=pow(e_at,0.25);
	B_lw=1+pow(e_at,0.25)+pow(e_at,0.5)+pow(e_at,0.75);
	/*Init_constant*/
	if (indicator_city==1)
	{
		Co_soil=Deep_land_c*CP_land_c;
		Co_solar_sw1=(1-Alpha_land_c)*(1-C_A)*(1-Alpha_atm)*Solar_const;
		Co_solar_sw_sin=Co_solar_sw1*sin(angle_lati)*sin(angle_dec);
		Co_solar_sw_cos=Co_solar_sw1*cos(angle_lati)*cos(angle_dec);
		Co_lw=Sky_VF_c*Sky_emi_c*Blotz*pow(ave_temp_2m,3.0)*B_lw;
	}
	else
	{
		Co_soil=Deep_land_r*CP_land_r;
		Co_solar_sw1=(1-Alpha_land_r)*(1-C_A)*(1-Alpha_atm)*Solar_const;
		Co_solar_sw_sin=Co_solar_sw1*sin(angle_lati)*sin(angle_dec);
		Co_solar_sw_cos=Co_solar_sw1*cos(angle_lati)*cos(angle_dec);
		Co_lw=Sky_VF_c*Sky_emi_r*Blotz*pow(ave_temp_2m,3.0)*B_lw;
	}
	
	Co_sh=9.8+10.5*U_ref;
	Co_up=Co_lw+Co_sh;
	Co_cons_1=-1*(Co_lw*B_e+Co_sh)*ave_temp_2m;
	Co_cons_2=Co_cons_1-Co_solar_sw_sin;
	SW=Co_solar_sw_sin-Co_solar_sw_cos*cos(omiga*time);

	if (SW<0.00001 && (time-time_start)<12*secs && start_hour*secs<=time_rise)
	{
		tmp_exp=exp(-1*Co_up*(time-time_start)/(Co_soil*secs));
		surf_temp=(Temp_land_ini+Co_cons_1/Co_up)*tmp_exp-Co_cons_1/Co_up;
		/*Message("Phase 1 at time %g is %g\n",time,surf_temp);*/
	}
	else if (SW>0.00001)
	{
		t_r=time_rise;
		tmp_exp_1=exp(-1*Co_up*time_rise/(Co_soil*secs));
		temp_ref_1=(Temp_land_ini+Co_cons_1/Co_up)*tmp_exp_1-Co_cons_1/Co_up;
		tmp_exp=exp(-1*(time-t_r)*Co_up/(Co_soil*secs));
		tmp_r=(Co_solar_sw_cos*Co_up*cos(omiga*t_r)+Co_soil*Co_solar_sw_cos*omiga*3600*sin(omiga*t_r))/(pow(Co_soil*omiga*3600,2.0)+Co_up*Co_up);
		tmp_n=(Co_solar_sw_cos*Co_up*cos(omiga*time)+Co_soil*Co_solar_sw_cos*omiga*3600*sin(omiga*time))/(pow(Co_soil*omiga*3600,2.0)+Co_up*Co_up);
		surf_temp=tmp_exp*(temp_ref_1+tmp_r+Co_cons_2/Co_up)-Co_cons_2/Co_up-tmp_n;
		/*Message("Phase 2 at time %g is %g\n",time,surf_temp);*/
	}
	else 
	{			
		t_r=time_set;
		tmp_exp_1=exp(-1*Co_up*time_rise/(Co_soil*secs));
		temp_ref_1=(Temp_land_ini+Co_cons_1/Co_up)*tmp_exp_1-Co_cons_1/Co_up;
		tmp_exp_2=exp(-1*(time_set-time_rise)*Co_up/(Co_soil*secs));
		tmp_r=(Co_solar_sw_cos*Co_up*cos(omiga*time_rise)+Co_soil*Co_solar_sw_cos*omiga*3600*sin(omiga*time_rise))/(pow(Co_soil*omiga*3600,2.0)+Co_up*Co_up);
		tmp_n=(Co_solar_sw_cos*Co_up*cos(omiga*time_set)+Co_soil*Co_solar_sw_cos*omiga*3600*sin(omiga*time_set))/(pow(Co_soil*omiga*3600,2.0)+Co_up*Co_up);
		temp_ref_2=tmp_exp_2*(temp_ref_1+tmp_r+Co_cons_2/Co_up)-Co_cons_2/Co_up-tmp_n;
		tmp_exp=exp(-1*Co_up*(time-t_r)/(Co_soil*secs));
		surf_temp=(temp_ref_2+Co_cons_1/Co_up)*tmp_exp-Co_cons_1/Co_up;
		/*Message("Phase 3 at time %g is %g\n",time,surf_temp);*/
	}	
	/*Message("surface temperature is : %f (K)\n",surf_temp);*/
#endif /*RP_HOST*/	
	return surf_temp;
}

/*******************************************************/
/*UDF functions*/
/********  UDM define *********/
DEFINE_ON_DEMAND(Init_UDM) /* run this first to initialize your UDM*/
{ 
/* Variables used by serial, host, node versions */
	real y,x_s,tmp_por;
	real tmp_K1,tmp_K2,tmp_C1,tmp_C2,tmp_beta;
	real Q_anth_c;	
/* "Parallelized" Sections */
#if !RP_HOST /* Host will do nothing in this udf. Serial will */
	Domain *domain;   /* declare domain pointer since it is not passed a argument to DEFINE macro  */
	real  x[ND_ND];
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
					tmp_por = Por_0;
					C_UDMI(c,t,1) = tmp_por;
					tmp_beta=1.0;
				}
				else if (x_s>=Por_xs2 && x_s<=Por_xe2)
                                {
                                        tmp_por = Por_0;
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

	surf_temp_c=surface_temp(1,ave_temp_2m_c,U_ref_c);
        surf_temp_r=surface_temp(0,ave_temp_2m_r,U_ref_r);
        Q_anth_c=anthro_heat(start_hour*secs);
		Q_value_sen_build=(9.8+4.1*U_ref_c)*(surf_temp_c-ave_temp_2m_c);
        Q_value_sen_city=Q_value_sen_build+Q_anth_c;
        Q_value_sen_rural=(9.8+4.1*U_ref_r)*(surf_temp_r-ave_temp_2m_r);
        #endif /* !RP_NODE */
	node_to_host_real_2(Q_value_sen_city,Q_value_sen_rural); /* Does nothing in SERIAL */
		#if !RP_NODE  /* SERIAL or HOST */
			Message("\nsensible heat flux in city is %g and in rural is %g\n",Q_value_sen_city,Q_value_sen_rural);
		#endif /* !RP_NODE */
}
DEFINE_ON_DEMAND(Init_FUDM)
{ 
	int surface_thread_id=9;
	real Q_anth_c; 	
#if !RP_HOST
	Domain *domain;   /* declare domain pointer since it is not passed a argument to DEFINE macro  */
	cell_t c0;
	Thread *t,*t0;
	face_t f;
#endif /* !RP_HOST */
	host_to_node_int_1(surface_thread_id);
#if !RP_HOST
	domain=Get_Domain(1); 
	t = Lookup_Thread(domain,surface_thread_id);
	begin_f_loop(f,t)
	{
		if (PRINCIPAL_FACE_P(f,t))
		{
			t0=THREAD_T0(t);
			c0=F_C0(f,t);
			F_UDMI(f,t,0)= C_T(c0,t0);
		}
	}
	end_f_loop(f,t)
#endif /* !RP_HOST */
}
/*calculate the surface temperature*/
DEFINE_ON_DEMAND(test_anthro)
{
	/*real time,time_start;
	time_start=start_hour*secs;
	time = CURRENT_TIME*iter_factor+time_start;*/
	/*Message("\ncall anthropogetic heat function at time %g, the heat is %g",time, anthro_heat(time));*/
	ave_temp_2m_c=(average_temp_zone(line_2m,eps_2m,Por_xs,Por_xe)+average_temp_zone(line_2m,eps_2m,Por_xs2,Por_xe2))/2.0;
   ave_temp_2m_r=average_temp_zone(line_2m,eps_2m,rural_r_s,rural_r_e);
   U_ref_c=(average_abs_vel_zone(line_10m,eps_10m,Por_xs,Por_xe)+average_abs_vel_zone(line_10m,eps_10m,Por_xs2,Por_xe2))/2.0;
   U_ref_r=average_abs_vel_zone(line_10m,eps_10m,rural_r_s,rural_r_e);
   Message("\n city average temperature %g, average velocity %g ",ave_temp_2m_c, U_ref_c);
   Message("\n rural average temperature %g, average velocity %g ",ave_temp_2m_r, U_ref_r);
}
DEFINE_EXECUTE_AT_END(Surf_temp)
{
	real time,time_start,Q_anth_c;
	time_start=start_hour*secs;
	time = CURRENT_TIME*iter_factor+time_start;	
	if (time/updata_time_tl==ceil(time/updata_time_tl))
	{	
		ave_temp_2m_c=(average_temp_zone(line_2m,eps_2m,Por_xs,Por_xe)+average_temp_zone(line_2m,eps_2m,Por_xs2,Por_xe2))/2.0;
		ave_temp_2m_r=average_temp_zone(line_2m,eps_2m,rural_r_s,rural_r_e);
		U_ref_c=(average_abs_vel_zone(line_10m,eps_10m,Por_xs,Por_xe)+average_abs_vel_zone(line_10m,eps_10m,Por_xs2,Por_xe2))/2.0;
		U_ref_r=average_abs_vel_zone(line_10m,eps_10m,rural_r_s,rural_r_e);
		surf_temp_c=surface_temp(1,ave_temp_2m_c,U_ref_c);
		surf_temp_r=surface_temp(0,ave_temp_2m_r,U_ref_r);
		Q_anth_c=anthro_heat(time);
		node_to_host_real_6(ave_temp_2m_c,U_ref_c,ave_temp_2m_r,U_ref_r,surf_temp_c,surf_temp_r); /* Does nothing in SERIAL */
	#if !RP_NODE  /* SERIAL or HOST */
		Message("\nAt time %g s, in city:\n",(CURRENT_TIME*iter_factor+time_start));
		Message("The air temperature at 2m is %g, the horizontal velicity at 10m is %g\n",ave_temp_2m_c, U_ref_c);
		Message("The surface temperature is %g\n",surf_temp_c);
		Message("In rural:\n");
		Message("The air temperature at 2m is %g, the horizontal velicity at 10m is %g\n",ave_temp_2m_r, U_ref_r);
		Message("The surface temperature is %g\n",surf_temp_r);
	#endif /* !RP_NODE */
		Q_value_sen_build=(9.8+4.1*U_ref_c)*(surf_temp_c-ave_temp_2m_c);
        Q_value_sen_city=Q_value_sen_build+Q_anth_c;
		Q_value_sen_rural=(9.8+4.1*U_ref_r)*(surf_temp_r-ave_temp_2m_r);
		node_to_host_real_2(Q_value_sen_city,Q_value_sen_rural);
	#if !RP_NODE  /* SERIAL or HOST */
			Message("HEATFLUX in city is %g and in rural is %g, for building is %g\n",Q_value_sen_city,Q_value_sen_rural,Q_value_sen_build);
	#endif /* !RP_NODE */
	}
}
DEFINE_EXECUTE_AT_END(temp_update)
{
	int surface_thread_id=9;
	real Q_anth_c; 	
#if !RP_HOST
	Domain *domain;   /* declare domain pointer since it is not passed a argument to DEFINE macro  */
	cell_t c0;
	Thread *t,*t0;
	face_t f;
#endif /* !RP_HOST */
	host_to_node_int_1(surface_thread_id);
#if !RP_HOST
	domain=Get_Domain(1); 
	t = Lookup_Thread(domain,surface_thread_id);
	begin_f_loop(f,t)
	{
		if (PRINCIPAL_FACE_P(f,t))
		{
			t0=THREAD_T0(t);
			c0=F_C0(f,t);
			F_UDMI(f,t,0)= C_T(c0,t0);
		}
	}
	end_f_loop(f,t)
#endif /* !RP_HOST */
}
/*****************************  profiles****************** *********/
/***************************porous************************/
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
/***************************heat flux*********************/
/******** city sensible heat flux  *********/
DEFINE_PROFILE(WallFlux_build,thread,i)
{	
	real source;
	face_t f;
	
	begin_f_loop(f,thread)
	{
		source=Q_value_sen_build*iter_factor;
		F_PROFILE(f,thread,i) = source;
	}
	end_f_loop(f,thread)

}
DEFINE_PROFILE(WallFlux_city,thread,i)
{	
	real source;
	face_t f;
	
	begin_f_loop(f,thread)
	{
		source=Q_value_sen_city*Por_0*iter_factor;
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
		source=Q_value_sen_rural*iter_factor;
		F_PROFILE(f,thread,i) = source;
	}
	end_f_loop(f,thread)

}
DEFINE_PROFILE(Temp_p_out,thread,i)
{	
	real source;
	face_t f;
	cell_t c0;
	Thread *t0;
	begin_f_loop(f,thread)
	{
		/*source=ave_temp_2m_r;*/
		F_PROFILE(f,thread,i) = F_UDMI(f,thread,0);
		/*Message("\n temporaly temperature is %g\n",F_UDMI(f,thread,0));*/
	}
	end_f_loop(f,thread)

}
/********  source term *********/
DEFINE_SOURCE(u_damping,c,t,dS,eqn)
{
	#if !RP_HOST
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
		damping_factor=-1*(1-cos(z_ref_height*pi))/2.0;
	}
	else 
	{
		damping_factor=-1*(1+(z_ref_height-0.5)*pi)/2.0;
	}
	con =damping_factor*C_U(c,t);
	source = con;
	dS[eqn] = damping_factor;
	return source;
	#endif
}
DEFINE_SOURCE(v_damping,c,t,dS,eqn)
{
	#if !RP_HOST
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
		damping_factor=-1*(1-cos(z_ref_height*pi))/2.0;
	}
	else 
	{
		damping_factor=-1*(1+(z_ref_height-0.5)*pi)/2.0;
	}
	con =damping_factor*C_V(c,t);
	J=C_UDMI(c,t,0);
	con2=Den_air*(J*J-1)*C3*g_acc*(C_T(c,t)-Temp_air_init);
	source = con+con2;
	dS[eqn] = damping_factor;
	return source;
	#endif
}
DEFINE_SOURCE(energy_source,c,t,dS,eqn)
{
	#if !RP_HOST
	real x[ND_ND];
	real y,con, source;
	real damping_factor;
	real z_ref_height;
	real con_por;
	real J;
	real Q_sen_city;
	C_CENTROID(x,c,t);
	y=x[1];
	z_ref_height=(y-Z_damping)/(Z_top-Z_damping);
	if (z_ref_height<=0)
	{
		damping_factor=0;
	}
	else if (z_ref_height<=0.5 && z_ref_height>0)
	{
		damping_factor=-1*(1-cos(z_ref_height*pi))/2.0;
	}
	else 
	{
		damping_factor=-1*(1+(z_ref_height-0.5)*pi)/2.0;
	}
	J=C_UDMI(c,t,0);
	con = C2*Den_air*CP_air*C_V(c,t)*J;
	if (Por_0==1.0)
	{
		source = con+damping_factor*(C_T(c,t)-Temp_air_init);
	}
	else
	{
		con_por=(1-Por_0)*Q_value_sen_city/H_Urb;
		source = con+damping_factor*(C_T(c,t)-Temp_air_init)+con_por;
	}
	dS[eqn] = damping_factor;
	return source;
	#endif
}
/****		TKE   source				****/
DEFINE_SOURCE(kinetic_source,c,t,dS,eqn)
{
	#if !RP_HOST
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
	if (Por_0==1.0)
	{
		source = con_stra;
		dS[eqn] = 0;
	}
	else
	{
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
	}
	return source;
	#endif
}
/****		TDR   source				****/
DEFINE_SOURCE(disspation_source,c,t,dS,eqn)
{
	#if !RP_HOST
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
	if (Por_0==1.0)
	{
		source = con_stra*ed;
		dS[eqn] = con_stra;
	}
	else
	{
		tmp_K=C_UDMI(c,t,1);
		tmp_C2=C_UDMI(c,t,2);;
		Q=sqrt(u*u+v*v);
		con1=-2*Por_0*mu_l*tmp_K;
		con2= (-4/3)*pow(Por_0,2)*tmp_C2*rho*Q;
		con=con1*ed+con2*ed;
		source = con+con_stra*C_D(c,t);
		dS[eqn] = con1+con2+con_stra;
	}
	return source;
	#endif
}
