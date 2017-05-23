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

#define MODELSRC "muPillar_FB.cc"
#define MODELPATH "muPillar_FB/"
#define PI 3.1415


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
	std::ofstream out;
	if(!katana::getCmdOption_bool(argv, argv+argc, "-nolog", false))
		{if(katana::getCmdOption_bool(argv, argv+argc, "-rmlog", false)) katana::rm(path+"log.txt");
		out.open(path+"log.txt", std::ios::app);
		std::cout.rdbuf(out.rdbuf());
		};
	//--- END Change cout to logfile --- //
	
	
	std::string sweepfilename=katana::getCmdOption(argv, argv+argc, "-sfn", "sweep.dat");
	bool append_data 	=katana::getCmdOption_bool(argv, argv+argc, "-append", false);
	std::ios::openmode iohandler; if(append_data) iohandler=std::ios::app; else iohandler=std::ios::trunc; 


	std::string tsfilename=katana::getCmdOption(argv, argv+argc, "-tsname", "ts.data");
	std::string turnontsfilename=katana::getCmdOption(argv, argv+argc, "-turnontsname", "turnonts.data");
	std::string g2filename=katana::getCmdOption(argv, argv+argc, "-g2filename", "g2.data");
	std::string specfilename=katana::getCmdOption(argv, argv+argc, "-specfilename", "fourier.data");
	std::string sourcecode	=katana::getCmdOption(argv, argv+argc, "-srccode", MODELSRC);


	//---- END PROGRAM DEPENDEND

	//----- CP Source Code Main to data/date Folder --- //
	katana::cp_src(sourcecode, path);
	//----- END CP ------------------------------------ //

//-----	END CONFIGURATION ---------- //

//----- BEGIN PARAMETERS

	Time t_start = 	katana::getCmdOption(argv, argv+argc, "-t_start", 0.)*1e-9;
	Time t_stop=	katana::getCmdOption(argv, argv+argc, "-t_stop", 50)*1e-9;
	Time t_turnon=	katana::getCmdOption(argv, argv+argc, "-t_turnon", 10)*1e-9;
	Time t_trysteady=katana::getCmdOption(argv, argv+argc, "-t_trysteady", 15)*1e-9;
	Time dt = 	katana::getCmdOption(argv, argv+argc, "-dt", 0.001)*1e-9;// calc in s, measure ns
	
	double J=0;	
	double jstart=katana::getCmdOption(argv, argv+argc, "-jstart" , 60.0)*1e-6; 			
	double jstop=katana::getCmdOption(argv, argv+argc, "-jstop" , 300)*1e-6; 			
	double jincr=katana::getCmdOption(argv, argv+argc, "-jincr" , 10.0)*1e-6; //Input in numbers of the Jth

	size_t corrprlen=int(katana::getCmdOption(argv, argv+argc, "-corrprintlen", 0)/dt);	
	size_t fourierlen=int(katana::getCmdOption(argv, argv+argc, "-fourierprintlen", 10));	
	size_t cutfouriertrace=int(katana::getCmdOption(argv, argv+argc, "-cutfouriertrace", 20));
	int makeavgs=katana::getCmdOption(argv, argv+argc, "-makeavgs", 1);

//---- END PARAMETERS

//----- VARS & INIT CONDITIONS ----------- // 
	//-- VARS 
 
	State X(6, 0);
	Time tau_fb = katana::getCmdOption(argv, argv+argc, "-tau_fb", 0.8)*1e-9;
	double	historyomega=katana::getCmdOption(argv, argv+argc, "-historyomega", 0);
	
	//-- INIT COND
	
	X[0]=katana::getCmdOption(argv,argv+argc, "-Esinit_Re", 0.1);
	X[1]=katana::getCmdOption(argv,argv+argc, "-Esinit_Im", 0.1);
	X[2]=katana::getCmdOption(argv,argv+argc, "-Ewinit_Re", 0.1);
	X[3]=katana::getCmdOption(argv,argv+argc, "-Ewinit_Im", 0.1);
	X[4]=katana::getCmdOption(argv,argv+argc, "-Rhoinit", 0.1);
	X[5]=katana::getCmdOption(argv,argv+argc, "-Ninit", 0.1);

	const size_t hist_len = int(tau_fb/dt)==0 ?  1 : int(tau_fb/dt); // HIER unnützer Overhead für Fall das ohne Delay integriert wird
	Delayvector X_tau_vector(hist_len, State (4,0));  //size MUST be initiliazed and equal (dt*tau) = indirect definition of tau 

//----- END VARS & INIT CONDITIONS ------ //


//----- MODEL PARAMS --------------------//
	//-- define params in kg, m, s
	//-- general constants, these are the last parameters passed to the par vector
	auto epsi0 = 8.85e-12;		//F m^-1 == J V^-2 m^-1
	auto hbar = 1.0545e-34;		//Js
	auto e0 = 1.6e-19;			//C
	auto c0 = 3.0e8;			//m s^-1
	std::complex<double>	ii={0,1};
	//-- Laser parameters
	auto kappa_s    = katana::getCmdOption(argv,argv+argc, "-kappa_s", 39.)*(1/1e-9);	//ps^-1 -> 1/s
	auto kappa_w    = katana::getCmdOption(argv,argv+argc, "-kappa_w", 42.)*(1/1e-9);	//ps^-1 -> 1/s
	auto mu_s 	= 3.70*(1e-9*e0);	//nm * e0 -> m*C
	auto mu_w 	= 3.75*(1e-9*e0);	//nm * e0 -> m*C
	auto _epsi_ss   = katana::getCmdOption(argv, argv+argc, "-gss", 70.)*1e-10;		//m^2 A^-1 V-^1
	auto _epsi_ww   = katana::getCmdOption(argv, argv+argc, "-gww", 50.)*1e-10;		//m^2 A^-1 V^-1
	auto _epsi_sw   = katana::getCmdOption(argv, argv+argc, "-gsw", 160.)*1e-10;		//m^2 A^-1 V^-1
	auto _epsi_ws   = katana::getCmdOption(argv, argv+argc, "-gws", 150.)*1e-10;		//m^2 A^-1 V^-1
	auto beta	= katana::getCmdOption(argv, argv+argc, "-beta", 5.6e-3);
	auto eta 	= katana::getCmdOption(argv, argv+argc, "-eta", 1.28)*1e-3;
	auto tau_r 	= katana::getCmdOption(argv, argv+argc, "-tau_r", 150.)*(1e-12); 	//ps -> s
	auto S_in 	= katana::getCmdOption(argv, argv+argc, "-S_in", 1.0e-4);//m^2 ps^-1 -> m^2 s^-1
	auto V	 	= 6.3*1E-18;// 0.5*6.3*(1E-6)*1E-6*1E-6; 	//micro m^3 -> m^3
	auto A		= 3.14* 1E-12; // 3.14*(1E-6)*1E-6; 	//micro m^2
	auto Z_QD 	= katana::getCmdOption(argv, argv+argc, "-ZQD", 110.);
	auto Jp		= katana::getCmdOption(argv, argv+argc, "-Jp", 42.5)*(1e-6);
	auto n_bg 	= 3.34;
	auto epsi_bg	= n_bg*n_bg;
	auto tau_sp 	= katana::getCmdOption(argv, argv+argc, "-tau_sp", 1.0)*(1e-9);		//ns -> s
	auto T_2 	= 0.33*(1e-12); 	//ps -> s
	auto hbar_omega	= 1.38*(1.6e-19); 	//eV -> J
	auto _epsi_tilda = epsi0*n_bg*c0;
	auto alpha	=0.;

	// Abbreviations
	auto PF_E=0.5*(hbar_omega/(epsi0*epsi_bg))*2.0*Z_QD/V;
	auto eps_ss =_epsi_tilda*_epsi_ss; 	auto eps_sw =_epsi_tilda*_epsi_sw;
	auto eps_ws =_epsi_tilda*_epsi_ws;	auto eps_ww =_epsi_tilda*_epsi_ww;
	auto pump_factor=eta/(e0*A);

	//continued param in paper 	//feedback/delay params
	auto feed_phase_ss = katana::getCmdOption(argv, argv+argc, "-fP_ss", 0.); 	auto feed_ampli_ss = katana::getCmdOption(argv, argv+argc, "-fA_ss", 0.);
	auto feed_phase_ww = katana::getCmdOption(argv, argv+argc, "-fP_ww", 0.);	auto feed_ampli_ww = katana::getCmdOption(argv, argv+argc, "-fA_ww", 0.);
	auto f_ss=feed_ampli_ss*std::exp(ii*feed_phase_ss); 				auto f_ww=feed_ampli_ww*std::exp(ii*feed_phase_ww);	 			//complex feedback strength
	double deltaOmega_s=0, deltaOmega_w=0;
	double deltaOmega_ACC=katana::getCmdOption(argv, argv+argc, "-deltaOm_ACC", 1E-3);				//Berechne den Rotating frame bis auf diese relative Genauigkeit
	bool sponE=katana::getCmdOption_bool(argv, argv+argc, "-sponE", false);		// Toggle spontaneous emission
	katana::Seed s;		s.set_seed();

//----- END MODEL PARMS -----------------//

//----- BEGIN FILEHEADERS ---------------//
	
	std::stringstream traceheader_stringstream;
	traceheader_stringstream << "#Time [" << tau_sp*1e9 << "ns] \t Re(Es) \t Im(Es) \t |Es|^2 \t rho \t w";
	std::string traceheader(traceheader_stringstream.str());
	std::string sweepfileheader="\n#J \t |Es|^2 \t |Ew|^2 \t rho \t w ";
	
//----- END FILEHEADERS -----------------//


//----- BEGIN FILL DELAYVECTOR IF NEEDED //
	if(katana::getCmdOption_bool(argv, argv+argc, "-fillhistory", false) || katana::getCmdOption_bool(argv, argv+argc, "-historyomega", false))
	 {
	  for(uint iter=0; iter<hist_len; ++iter)
		{ // Retrieve Amp and Phase from X[0]:
			double Amp_s=std::sqrt(X[0]*X[0]+X[1]*X[1]);
			double Phase_s=std::arg(std::complex<double> {X[0],X[1]});
			double Amp_w=std::sqrt(X[2]*X[2]+X[3]*X[3]);
			double Phase_w=std::arg(std::complex<double> {X[2],X[3]});
			X_tau_vector[iter][0]=Amp_s*std::cos(historyomega/(2*PI)*iter*dt/tau_sp-Phase_s);
			X_tau_vector[iter][1]=Amp_s*std::sin(historyomega/(2*PI)*iter*dt/tau_sp-Phase_s);
			X_tau_vector[iter][2]=Amp_w*std::cos(historyomega/(2*PI)*iter*dt/tau_sp-Phase_w);
			X_tau_vector[iter][3]=Amp_w*std::sin(historyomega/(2*PI)*iter*dt/tau_sp-Phase_w);
			X_tau_vector[iter][4]=X[4];
			X_tau_vector[iter][5]=X[5];
		}
          }
	 else std::fill(X_tau_vector.begin(), X_tau_vector.end(), X);
//----- END FILL DELAYVECTOR ------------//


//----- BEGIN LAMBDA FUNCTION FOR DIFFERENTIAL EQUATION, MAPS X to dxdt ---- //
	auto dde_equation = [&](State &X, State &X_tau, Time &T, State &dxdt) 
		{
			// disassembling State vector into something nice
			
			std::complex<double>	Es={X[0],X[1]},	Ew={X[2],X[3]};
			std::complex<double>	Es_tau={X_tau[0],X_tau[1]},Ew_tau={X_tau[2],X_tau[3]};
			auto rho=X[4], w=X[5];

			std::complex<double>	dxdtEs={0,0}, dxdtEw={0,0};
			std::complex<double>	dxdtEs_sp, dxdtEw_sp;
			double dxdtrho=0, dxdtw=0;
			auto	gwN1=katana::gwNoise(), gwN2=katana::gwNoise();

			if(sponE) {dxdtEs_sp=std::sqrt(beta*rho/tau_sp*PF_E/dt)*gwN1;} else dxdtEs_sp=beta*rho/tau_sp*PF_E/std::conj(Es);
			if(sponE) {dxdtEw_sp=std::sqrt(beta*rho/tau_sp*PF_E/dt)*gwN2;} else dxdtEw_sp=beta*rho/tau_sp*PF_E/std::conj(Ew);
			auto gs = (std::norm(mu_s/hbar))*T_2*0.5/(1.+eps_ss*std::norm(Es)+ eps_sw*std::norm(Ew));
			auto gw = (std::norm(mu_w/hbar))*T_2*0.5/(1.+eps_ws*std::norm(Es)+eps_ww*std::norm(Ew));
			//double rhoinac=tau_sp*S_in*w/(1+tau_sp*S_in*w); //tau_sp*tau_sp*R*R*w/(1+tau_sp*tau_sp*R*R*w);
			double S=S_in;
			

			dxdtEs=PF_E*gs*(2*rho-1)*(double(1)-ii*alpha)*Es-kappa_s*(Es-f_ss*Es_tau)+dxdtEs_sp-ii*deltaOmega_s*Es;
			dxdtEw=PF_E*gw*(2*rho-1)*(double(1)-ii*alpha)*Ew-kappa_w*(Ew-f_ww*Ew_tau)+dxdtEw_sp-ii*deltaOmega_w*Ew;
			dxdtrho=-(2*rho-1)*(gs*std::norm(Es)+gw*std::norm(Ew))+S*w*(1-rho)-rho/tau_sp;
			dxdtw=pump_factor*(J-Jp)-S*2*Z_QD/A*w*(1-rho)- w/tau_r ;

			dxdt[0]=dxdtEs.real(); 	dxdt[1]=dxdtEs.imag();
			dxdt[2]=dxdtEw.real();	dxdt[3]=dxdtEw.imag();
			dxdt[4]=dxdtrho; 	dxdt[5]=dxdtw;
	
		};
//----- END LAMBDA FOR DIFFERENTIAL EQUATION ------------------------------ //	
	
	
//----- PRINT HEADER OF CONFIG/SWEEPFILE ---------------------------------- //
	sweepfilename=path+sweepfilename;	
	std::ofstream outputfile (sweepfilename, std::ios::out | iohandler);      //Output File (ofstream) overwrite
	
	outputfile << "# Model from RED16, control program to compare with BKM. project started: 2016/10/26 : last updated: 2016/10/26" << std::endl;
	outputfile << "# Parameters: to be found in source code and log.txt!" << std::endl; 
	outputfile << "#Columns:" << sweepfileheader  << std::endl;
	
//------ BEGIN PREPARATION ----------------------------//

	size_t simulation_steps=int((t_stop-t_start)/dt);
	int present_simustep=0;
	std::vector<std::complex<double>> Es_trace, Ew_trace;
	std::vector<double> Ints_trace, Intw_trace;
	std::vector<State> Full_Ttrace;
	katana::container_1d avgstate, avgacs, avgacw, avgcc, avgspecs, avgspecw, deltaOmcontainer, avg_Ints, avg_Intw;
	
//----- Lambda function for output, easier to implement inside the main() because we already have declared State, Time and Delayvector type		

	
	auto output_FullTtrace = [&](State &X, State &X_tau, Time &T, State &dxdt) 
		{	Full_Ttrace.push_back(X);
		};
		
	auto output_Int_Fields = [&](State &X, State &X_tau, Time &T, State &dxdt) 
		{	double Ints=X[0]*X[0]+X[1]*X[1], Intw=X[2]*X[2]+X[3]*X[3], rho=X[4], w=X[5];
			Ints_trace.push_back(Ints); Intw_trace.push_back(Intw);
			if(present_simustep%cutfouriertrace==0) { Es_trace.push_back({X[0], X[1]}); Ew_trace.push_back({X[2], X[3]}); present_simustep=0;}
			++present_simustep;
			avgstate.add({Ints, Intw, rho, w});
		};

	auto eval_deltaOm = [&](State &X, State &X_tau, Time &T, State &dxdt)
		{ 	std::complex<double> ImEEs=(dxdt[0]+ ii* dxdt[1])/(X[0]+ ii*X[1]);
			std::complex<double> ImEEw=(dxdt[2]+ ii* dxdt[3])/(X[2]+ ii*X[3]);
			std::vector<double> deltaOm={ImEEs.imag()+deltaOmega_s, ImEEw.imag()+deltaOmega_w};
			deltaOmcontainer.add(deltaOm);
		};

//---- DOING J-SWEEP  ------- //
	
	for(J=jstart; J<=jstop; J=J+jincr)
		{//---- BEGIN PREPARATION ------------------------------//
			std::string trace_outname, acs_outname, acw_outname, cc_outname, specs_outname, specw_outname, Intwtrace_outname, Intstrace_outname;
			char J_buf[10]; 	std::sprintf(J_buf , "%06.2f", J*1E6);
			trace_outname= path +  J_buf + "_muA_" + turnontsfilename; 
			Intstrace_outname= path +  J_buf + "_muA_s_" + tsfilename; 
			Intwtrace_outname= path +  J_buf + "_muA_w_" + tsfilename; 
			acs_outname= path +  J_buf + "_muA_" + "auto_s_" +g2filename;
			acw_outname= path +  J_buf + "_muA_" + "auto_w_" +g2filename;
			cc_outname= path +  J_buf + "_muA_" + "cross_sw_" +g2filename;
			specs_outname= path +  J_buf + "_muA_" + "optspec_s_" +specfilename;
			specw_outname= path +  J_buf + "_muA_" + "optspec_w_" +specfilename;
		
			std::cout << "Change current to J= "<< J*1E6 << std::endl;
			avgstate.clear(); avgacs.clear(); avgacw.clear(); avgcc.clear(), deltaOmcontainer.clear(), avg_Ints.clear(), avg_Intw.clear(); // Clear all container
			Full_Ttrace.clear(); 				// Making sure the array is clear
		//----- END PREPAREATION _------------------------------//
		//----- MEASURE TURNON ---------------------------------//
			roehm::ddesolve_euler_singledelay(dde_equation, X, X_tau_vector, t_start, t_turnon , dt, output_FullTtrace); // 5 ns to simulate turnon
			katana::dump_trace(Full_Ttrace,dt, trace_outname, traceheader, 10);
		//----- END TURNON ------------------------------------//

		//----  BEGIN Iterating to "steady state", calculating rotating frame within realtive accuarcy of deltaOmega_ACC -----//
			int trysteady=0;
			while(true)
				{ roehm::ddesolve_euler_singledelay(dde_equation, X, X_tau_vector, t_start, t_trysteady , dt, eval_deltaOm); // 5 ns to simulate turnon
				  std::vector<double> deltaOm_vector; deltaOmcontainer.get_avg(deltaOm_vector);
				  if(deltaOmega_ACC> std::abs(deltaOmega_s/deltaOm_vector[0]-1) && deltaOmega_ACC > std::abs(deltaOmega_w/deltaOm_vector[1]-1))
					{std::cout << "Reached steady state! Measure!" << std::endl;
					 std::cout << "Change rotating frame to "<< deltaOmega_s << " (Strong), " << deltaOmega_w << "(weak)" <<std::endl;
					 break;}
				  deltaOmega_s=deltaOm_vector[0]; deltaOmega_w=deltaOm_vector[1];
				  ++trysteady;
				  std::cout << "Did "<< trysteady << " Iterations by now. Not yet at steadystate, iterating once more!" << std::endl;
				  if(trysteady>100) { std::cout<< "Enough trying, take average rotating frame after 100 iterations and hope the best" << std::endl; break;}
				}
		//---- END Iterate steady state --//

		//----- BEGIN MEASURING! ----------//
			for(int iter=0; iter<makeavgs; ++iter)
				{std::cout << "Start Iteration: " << iter+1<< " of " << makeavgs << std::endl;
				Ints_trace.clear();	Intw_trace.clear(); 	// Making sure the array is clear
				Es_trace.clear();	Ew_trace.clear(); 	// Making sure the array is clear
				Full_Ttrace.reserve(simulation_steps);
				Ints_trace.reserve(simulation_steps); 	Intw_trace.reserve(simulation_steps);
				Es_trace.reserve(simulation_steps);	Ew_trace.reserve(simulation_steps);
				roehm::ddesolve_euler_singledelay(dde_equation, X, X_tau_vector, t_start, t_stop, dt, output_Int_Fields);
				std::vector<double> opticalspec_s, opticalspec_w;
				katana::get_optical_spec(Es_trace, opticalspec_s); avgspecs.add(opticalspec_s);
				katana::get_optical_spec(Ew_trace, opticalspec_w); avgspecw.add(opticalspec_w);
				std::vector<double> acs, acw, cc; 
				katana::get_g2(Ints_trace, acs);		avgacs.add(acs);
				katana::get_g2(Intw_trace, acw);		avgacw.add(acw);
				katana::get_g2(Ints_trace, Intw_trace, cc);	avgcc.add(cc);
				// avg_Ints.add(Ints_trace); avg_Intw.add(Intw_trace);
				};							
		//----- END MEASURING ----//

		//----- BEGIN OUTPUT -----//
			std::vector<double> avg_vec;  avgstate.get_avg(avg_vec); 		std::vector<double> avg_cc;  avgcc.get_avg(avg_cc);
			std::vector<double> avg_acs;  avgacs.get_avg(avg_acs); 			std::vector<double> avg_acw;  avgacw.get_avg(avg_acw);
			std::vector<double> avg_outspec_s; avgspecs.get_avg(avg_outspec_s);	std::vector<double> avg_outspec_w; avgspecw.get_avg(avg_outspec_w);
			//std::vector<double> avg_Intw_vec; avg_Intw.get(avg_Intw_vec); std::vector<double> avg_Ints_vec; avg_Ints.get(avg_Ints_vec); 

			outputfile << J*1E6 << '\t' << avg_vec[0] << '\t' << avg_vec[1] << '\t' << avg_vec[2] << '\t' << avg_vec[3] << std::endl;
			katana::print_correlation(avg_acs ,dt*1e9, acs_outname, corrprlen, 5, "ns");
			katana::print_correlation(avg_acw ,dt*1e9, acw_outname, corrprlen, 5, "ns");
			katana::print_correlation(avg_cc ,dt*1e9, cc_outname, corrprlen, 5, "ns");
			katana::dump_trace(Ints_trace,dt, Intstrace_outname, "#Time[ns], Ints [S.I]", 10);
			katana::dump_trace(Intw_trace,dt, Intwtrace_outname, "#Time[ns], Intw [S.I]", 10);
			katana::dump_optical_spec(avg_outspec_s, cutfouriertrace*dt*1e9, (int) fourierlen*(t_stop-t_start)*1e9, specs_outname); //most of the spectrum is shit, only first X elem. dumped
			katana::dump_optical_spec(avg_outspec_w, cutfouriertrace*dt*1e9, (int) fourierlen*(t_stop-t_start)*1e9, specw_outname);
//---------------------END OUTPUT ------//
		}
//------END SWEEP! ---//	
	std::cout << "Done!" << std::endl;
	std::cout.rdbuf(cout_oss);
	return 0;
}


