#include <cstdlib>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <complex>
#include <climits>
#include <boost/filesystem.hpp>
#include "ddesolve_euler_singledelay.hh"
#include "inc_modules.hh"

#define MODELSRC "testben.cc"
#define MODELPATH "testben/"



int main(int argc, char** argv)
{	using State = std::vector<double>;
	using Time = double;
	using Delayvector = std::vector<State>; 
	
//------ BEGIN CONFIGURATION -------- //


	//---- BEGIN PROGRAM DEPENDEND
	std::string model	=MODELPATH;
	std::string path	=katana::getCmdOption(argv, argv+argc, "-path", "data/"); 		
	std::string date	=katana::getDateString(); 
 	std::string subpath	=katana::getCmdOption(argv, argv+argc, "-spath", "");
	path			=path+model+date+subpath;
	boost::filesystem::create_directories(path.c_str());

	//--- Change cout to logfile --- //
	std::streambuf *cout_oss=std::cout.rdbuf();	
	if(!katana::getCmdOption_bool(argv, argv+argc, "-nolog", false))
		{
		if(katana::getCmdOption_bool(argv, argv+argc, "-rmlog", false)) katana::rm(path+"log.txt");
		std::ofstream out(path+"log.txt", std::ios::app);
		std::cout.rdbuf(out.rdbuf());
		};
	//--- END Change cout to logfile --- //
	
	
	std::string sweepfilename=katana::getCmdOption(argv, argv+argc, "-sfn", "sweep.dat");
	bool append_data 	=katana::getCmdOption_bool(argv, argv+argc, "-append", false);
	std::ios::openmode iohandler; if(append_data) iohandler=std::ios::app; else iohandler=std::ios::trunc; 


	std::string tsfilename=katana::getCmdOption(argv, argv+argc, "-tsname", "_ts.data");
	std::string sourcecode	=katana::getCmdOption(argv, argv+argc, "-srccode", MODELSRC);

	//---- END PROGRAM DEPENDEND


	//----- CP Source Code Main to data/date Folder --- //
	katana::cp_src(sourcecode, path);
	//----- END CP ------------------------------------ //

//-----	END CONFIGURATION ---------- //

//----- BEGIN PARAMETERS

	double J=160E-6;
	
	Time t_start = 	katana::getCmdOption(argv, argv+argc, "-t_start", 0.);
	Time t_stop=	katana::getCmdOption(argv, argv+argc, "-t_stop", 50)*1E-9;
	Time dt = 	katana::getCmdOption(argv, argv+argc, "-dt", 0.1)*1E-12;// calc in s, measure ns

//---- END PARAMETERS

//----- VARS & INIT CONDITIONS ----------- // 
	//-- VARS 
 
	State X(4, 0);
	Time tau_fb = 6E-9;
	std::string sweepfileheader="\n#J \t ReE \t ImE \t rho \t w ";
	
	//-- INIT COND
	
	X[0]=0.1;
//	X[1]=.1;
	X[2]=0;
	X[3]=0;
	const size_t hist_len = (int)(tau_fb/dt)+1; // HIER unnützer Overhead für Fall das ohne Delay integriert wird
		
	Delayvector X_tau_vector(hist_len);  //size MUST be initiliazed and equal (dt*tau) = indirect definition of tau 

//----- END VARS & INIT CONDITIONS ------ //


//----- MODEL PARAMS --------------------//
	//-- define params in kg, m, s
	//-- general constants, these are the last parameters passed to the par vector
	auto epsi0 = 8.85e-12;		//F m^-1 == J V^-2 m^-1
	auto hbar = 6.58e-16;		//Js
	auto e0 = 1.6e-19;			//C
	auto c0 = 3.0e8;			//m s^-1
	//-- Laser parameters
	auto kappa_s    = 0.039*(1/1e-12);	//ps^-1 -> 1/s
//	auto kappa_w    = 0.041*(1/1e-12);	//ps^-1 -> 1/s
	auto mu_s 	   = 3.70*(1e-9*e0);	//nm * e0 -> m*C
//	auto mu_w 	   = 3.75*(1e-9*e0);	//nm * e0 -> m*C
	auto epsi_ss    = 70e-10;		//m^2 A^-1 V-^1
//	auto epsi_ww    = 50e-10;		//m^2 A^-1 V^-1
//	auto epsi_sw    = 160e-10;		//m^2 A^-1 V^-1
//	auto epsi_ws    = 150e-10;		//m^2 A^-1 V^-1
	auto beta	   = 5.6e-3;
	auto J_p 	   = 42.5*1e-6;		//microAmps -> Amps
	auto eta 	   = 1.28e-3;
	auto tau_r 	   = 150*(1e-12); 	//ps -> s
	auto S_in 	   = 1E-4;//m^2 ps^-1 -> m^2 s^-1
	auto V 	   = 6.3*(1E-6)*1E-6*1E-6; 	//micro m^3 -> m^3
	auto Z_QD 	   = 110;
	auto n_bg 	   = 3.34;
	auto epsi_bg	   = n_bg*n_bg;
	auto tau_sp 	   = 1*(1e-9);		//ns -> s
	auto T_2 	   = 0.33*(1e-12); 	//ps -> s
	auto A 	   = 3.14*(1E-6)*1E-6; 	//micro m^2
	auto hbar_omega = 1.38*(1.6e-19); 	//eV -> J
	auto epsi_tilda = epsi0*n_bg*c0;
	//continued param in paper
	//feedback/delay params
	auto feed_phase = 0.;
	auto feed_ampli = 0.05;


//----- END MODEL PARMS -----------------//

//----- BEGIN FILEHEADERS ---------------//
	
	std::stringstream traceheader_stringstream;
	traceheader_stringstream << "Time [" << tau_sp*1e9 << "ns] \t Re(Es) \t Im(Es) \t |Es|^2 \t rho \t w";
	std::string traceheader(traceheader_stringstream.str());
	
//----- END FILEHEADERS -----------------//


//----- BEGIN LAMBDA FUNCTION FOR DIFFERENTIAL EQUATION, MAPS X to dxdt ---- //
	auto dde_equation = [&](State &X, State &X_tau, Time &T, State &dxdt) 
		{
			// disassembling State vector into something nice
			
	std::complex<double>	Es=	 {X[0],X[1]};
	std::complex<double>	Es_tau= {X_tau[0], X_tau[1]};
	auto rho=X[2];
	auto w=X[3];

	auto complex_fb_strength_s= feed_ampli*(cos(feed_phase), -sin(feed_phase));  // calculate once, else significant speed loss!
	std::complex<double>	dxdtEs;
	double dxdtrho=0;
	double dxdtw=0;

//	auto gwN1=katana::gwNoise();
//	auto gwN2=katana::gwNoise();

	dxdtEs=(hbar_omega/(epsi0*epsi_bg))*(2*Z_QD/V)*(mu_s*mu_s*T_2/(2*hbar*hbar))* pow(1+epsi_ss*epsi_tilda*std::norm(Es), -1) *(2*rho-1)*Es-kappa_s*(Es-complex_fb_strength_s*Es_tau)+(beta*hbar_omega/(epsi0*epsi_bg))*(2*Z_QD/V)*(rho/tau_sp)*(Es/(std::norm(Es)));
	dxdtrho = -(std::norm(mu_s)*T_2/(2*hbar*hbar))*pow(1 + epsi_ss*epsi_tilda*std::norm(Es), -1)*(2*rho-1)*std::norm(Es)-rho/tau_sp+S_in*w*(1-rho);
	dxdtw = (eta/(e0*A))*(J-J_p)+S_in*w*(2*Z_QD/A)*(1-rho)-w/tau_r;    
	dxdt[0]=dxdtEs.real();
	dxdt[1]=dxdtEs.imag();
	dxdt[2]=dxdtrho;
	dxdt[3]=dxdtw;
	
		};
//----- END LAMBDA FOR DIFFERENTIAL EQUATION ------------------------------ //	
	
	
//----- PRINT HEADER OF CONFIG/SWEEPFILE ---------------------------------- //
	sweepfilename=path+sweepfilename;	
	std::ofstream outputfile (sweepfilename, std::ios::out | iohandler);      //Output File (ofstream) overwrite
	
	outputfile << "# Model from RED16, control program to compare with BKM. project started: 2016/10/26 : last updated: 2016/10/26" << std::endl;
	outputfile << "# Parameters: to find in source code and log.txt!" << std::endl; 
	outputfile << "#Columns:" << sweepfileheader  << std::endl;
	
//------ BEGIN PREPARATION ----------------------------//

	std::vector<std::complex<double>> E_trace;
	std::vector<State> Ttrace;
	
	double jstart=katana::getCmdOption(argv, argv+argc, "-jstart" , 60.0)*1E-6; 			
	double jstop=katana::getCmdOption(argv, argv+argc, "-jstop" , 200)*1E-6; 			
	double jincr=katana::getCmdOption(argv, argv+argc, "-jincr" , 10.0)*1E-6; //Input in numbers of the Jth

	katana::Seed s;
	s.set_seed();
	//Lambda function for output, easier to implement inside the main() because we already have declared State, Time and Delayvector type		
	auto output = [&](State &X, State &X_tau, Time &T, State &dxdt) 
		{
			Ttrace.push_back(X);
		};

//---- DOING J-SWEEP  ------- //
	
	for(J=jstart; J<=jstop; J=J+jincr)
		{

			//------ BEGIN PREPARATION ------------------------------//
			Ttrace.clear(); 	// Making sure the array is clear
			Ttrace.reserve((t_stop-t_start)/dt);
//			std::vector<double> outspec;
			roehm::ddesolve_euler_singledelay(dde_equation, X, X_tau_vector, t_start, t_stop, dt, output);
//			std::vector<double> opticalspec;
//			katana::get_optical_spec(E_trace, opticalspec);
			outputfile << J*1E6 << '\t' << X[0] << '\t' << X[1] << '\t' << X[2] << '\t' << X[3] << std::endl;
			std::string outname;
//
			char J_buf[10];
			std::sprintf(J_buf , "%.2f", J*1E6);
			outname= path +  J_buf + "_Jth_" + tsfilename;
			std::cout << outname << std::endl;
			katana::dump_trace(Ttrace,dt*1E9, traceheader, outname, 100);
//			katana::dump_optical_spec(opticalspec, dt, (int) fourierlength*(t_stop-t_start), outname); //most of the spectrum is shit, only first X elem. dumped
		}
	
	std::cout << "Done!" << std::endl;
	std::cout.rdbuf(cout_oss);
	return 0;
}


