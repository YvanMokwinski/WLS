#include "app.hpp"
#include "wls-io.hpp"

#include "wls-sparse-matrix.hpp"
#include "wls-direct-mkl-pardiso.hpp"

void wls_gesv_manpage(const char * filename_,
		      int * err_)
{
  err_[0] 	= 0;
  FILE* f 	= fopen(filename_,"w");
  if (!f)
    {
      
      err_[0] = 1;
      fprintf(stderr,"mokflow_make_manpage failed, can not open file '%s'\n",filename_);
      return;
    }
  fprintf(f,".TH wls-gesv \"1\" \"" __DATE__ "\"  \"wls-gesv description\"\n\n");
  fprintf(f,".SH NAME\n");
  fprintf(f,"\n");
  fprintf(f,"wls-gesv\n");
  fprintf(f,"\n");
  fprintf(f,".SH SYNOPSIS\n");
  fprintf(f,"\n");
  fprintf(f,".B wls-gesv <filename> [<options>]\n");
  fprintf(f,"\n");
  fprintf(f,".SH DESCRIPTION\n");
  fprintf(f,"\n");
  fprintf(f,".B wls-gesv \n");
  fprintf(f,"is a program to gesv the validity of a matrix described with the matrix market format.\n");
  fprintf(f,"\n");
  fprintf(f,".SH EXAMPLES\n");
  fprintf(f,"\n");
  fprintf(f,".SH AVAILABLE FILE FORMATS\n");
  fprintf(f,".BR\n");
  fprintf(f,".SS .mtx Matrix Market Format\n");
  fprintf(f,".br\n");
  fprintf(f,"see http://www\n");
  fprintf(f,"\n");
  fprintf(f,".SH OPTIONS\n");
  fprintf(f,"\n");
  fprintf(f,".SS\n");
  fprintf(f,".TP\n");
  fprintf(f,".B --help\n");
  fprintf(f,"Display help\n");
  fprintf(f,"\n");
  fprintf(f,".TP\n");
  fprintf(f,".B --version\n");
  fprintf(f,"Display the version number\n");
  fprintf(f,"\n");
  fprintf(f,".TP\n");
  fprintf(f,".B --manpage <filename>\n");
  fprintf(f,"Produce this manpage\n");
  fprintf(f,"\n");
  fprintf(f,".TP\n");
  fprintf(f,".B --nt <integer>\n");
  fprintf(f,"Set the number of threads\n");
  fprintf(f,"\n");
  fprintf(f,".TP\n");

#if 0
  fprintf(f,".B #######  STRING OPTIONS\n");
  fprintf(f,"\n");
  {  
    enum __emk_mesh_soption i = __emk_mesh_soption_error;
    for (++i;i<__emk_mesh_soption_n;++i)
      {
	const char * opt = __emk_mesh_soption_string(i);
	if (opt)
	  {
	    fprintf(f,".SS\n");
	    fprintf(f,".TP 5\n");
	    fprintf(f,".B %s <string>\n",opt);
	    fprintf(f,"%s\n",__emk_mesh_soption_desc(i));
	  }
      }
  }
  fprintf(f,".TP\n");
  fprintf(f,"\n\n.B #######  REAL OPTIONS\n");
  fprintf(f,"\n");
  {  
    enum __emk_mesh_roption i = __emk_mesh_roption_error;
    for (++i;i<__emk_mesh_roption_n;++i)
      {
	const char * opt = __emk_mesh_roption_string(i);
	if (opt)
	  {
	    fprintf(f,".SS\n");
	    fprintf(f,".TP 5\n");
	    fprintf(f,".B %s <real>\n",__emk_mesh_roption_string(i));
	    fprintf(f,"%s\n",__emk_mesh_roption_desc(i));
	  }
      }
  }
  fprintf(f,".TP\n");
  fprintf(f,"\n\n.B #######  INTEGER OPTIONS\n");
  fprintf(f,"\n");
  {  
    enum __emk_mesh_ioption i = __emk_mesh_ioption_error;
    for (++i;i<__emk_mesh_ioption_n;++i)
      {
	const char * opt = __emk_mesh_ioption_string(i);
	if (opt)
	  {
	    fprintf(f,".SS\n");
	    fprintf(f,".TP 5\n");
	    fprintf(f,".B %s <integer>\n",__emk_mesh_ioption_string(i));
	    fprintf(f,"%s\n",__emk_mesh_ioption_desc(i));
	  }
      }
  }  
  fprintf(f,".TP\n");
  fprintf(f,"\n\n.B #######  LOGICAL OPTIONS\n");
  fprintf(f,"\n");
  {  
    enum __emk_mesh_loption i = __emk_mesh_loption_error;
    for (++i;i<__emk_mesh_loption_n;++i)
      {
	const char * opt = __emk_mesh_loption_string(i);
	if (opt)
	  {
	    fprintf(f,".SS\n");
	    fprintf(f,".TP 5\n");
	    fprintf(f,".B %s\n",__emk_mesh_loption_string(i));
	    fprintf(f,"%s\n",__emk_mesh_loption_desc(i));
	  }
      }
  }  
#endif
  
  fprintf(f,"\n");
  fprintf(f,".TP\n");
  fprintf(f,".SH AUTHOR\n");
  fprintf(f,"Written by Yvan Mokwinski\n");
  err_[0] = fclose(f);

};



WLS::status_t calculate_solution(WLS::inverse_operator& 	inverse_op_,
				 const double * 		rhs_,
				 double * 		x_,
				 wls_int_t          rwork_n_,
				 double * 		rwork_)
{
  bool hasFailed;
  inverse_op_.compute(&hasFailed);
  if (hasFailed)
    {
      std::cerr << "compute failed: " << inverse_op_.get_error_message()  << std::endl;
      exit(1);
    }
  
  inverse_op_.apply("No transpose",
		    x_,
		    rhs_,
		    rwork_n_,
		    rwork_,
		    &hasFailed);
  
  if (hasFailed)
    {
      std::cerr << "apply failed: " << inverse_op_.get_error_message()  << std::endl;
      exit(1);
    }
  return WLS::status_t::success;
  
};

int main(int 		argc,
	 char ** 	argv)
{
  
  using status_t = WLS::status_t;
  
  WCOMMON::cmdline cmd(argc,argv);
  cmd.disp_header(stdout);
 
  if (cmd.get_logical("--help"))
    {
      fprintf(stdout,"// see manpage\n");
      return status_t::success;
    }

  {
    if (cmd.get_logical("--version"))
      {
	fprintf(stdout,"1.0\n");
	return status_t::success;
      }
  }

  {
    char filename[512];
    if (cmd.get_string("--manpage",filename))
      {
	int err;
	// "doc/man3/wls-gesv.3.gz"
	wls_gesv_manpage(filename, &err);
	return err;
      }
  }
  
  bool verbose = cmd.get_logical("-v");
  wls_int_t nt;
  if (!cmd.get_integer("--nt",&nt))
    {
      if (verbose)
	{
	  std::cerr << "// WARNING, missing '--nt <numthreads>', default is 1." << std::endl;
	}
      nt = 1;
    }
  
  char ofilename[512];
  if (!cmd.get_string("-o", ofilename))
    {
      std::cerr << "// WARNING, missing '-o <filename>', default is 'out.mtx'." << std::endl;
      sprintf(ofilename,"out.mtx");
    }
  
  if (cmd.get_nargs() < 3)
    {
      std::cerr << "ERROR missing files" << std::endl;
      return 1;
    }

  const char * a_filename = cmd.get_arg(1);
  const char * rhs_filename = cmd.get_arg(2);
  using real_t = double;
  
  WLS::sparse::matrix_t<real_t> a;
  status_t status = WLS::input::matrix_market_t::import(&a,
							a_filename);
  if (status)
    {
      return status;
    }

  
  WLS::dense::matrix rhs;
  status = WLS::input::matrix_market_t::import(rhs,
					       rhs_filename);
  if (status)
    {
      return status;
    }
  
  WLS::dense::matrix x(rhs.m,rhs.n);


  
  WLS::Direct::MKL::Pardiso pardiso(a.GetN(),a.GetB(),a.GetI(),a.GetX());

  double* rwork = nullptr;
  wls_int_t rwork_n = pardiso.get_buffer_size();
  if (rwork_n>0)
    {
      rwork = new double[rwork_n];
    }

  calculate_solution(pardiso,rhs(0,0),x(0,0),rwork_n,rwork);

  if (nullptr != rwork)
    {
      delete[] rwork;
      rwork = nullptr;
    }
  
  { WLS::output::matrix_market_t output(ofilename);
  output << "%\n";
  output << "% date: " << __DATE__ << "\n";  
  output << "% solution linear system, matrix: '" << a_filename << "', rhs: '" << rhs_filename << "'\n";
  output << "%\n";
  output << x; }
  //  std::cout << x << std::endl;
  return 0;
 
}

