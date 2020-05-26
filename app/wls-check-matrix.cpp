#include "wls.hpp"
#include "app.hpp"
#include "wls-io.hpp"

static void make_manpage(const char * filename_,
			      int * err_)
{
  err_[0] 	= 0;
  FILE* f 	= fopen(filename_,"w");
  if (!f)
    {
      err_[0] = 1;
      fprintf(stderr,"make_manpage failed, can not open file '%s'\n",filename_);
      return;
    }
  fprintf(f,".TH wls-check \"1\" \"" __DATE__ "\"  \"wls-check description\"\n\n");
  fprintf(f,".SH NAME\n");
  fprintf(f,"\n");
  fprintf(f,"wls-check\n");
  fprintf(f,"\n");
  fprintf(f,".SH SYNOPSIS\n");
  fprintf(f,"\n");
  fprintf(f,".B wls-check <filename> [<options>]\n");
  fprintf(f,"\n");
  fprintf(f,".SH DESCRIPTION\n");
  fprintf(f,"\n");
  fprintf(f,".B wls-check \n");
  fprintf(f,"is a program to check the validity of a matrix described with the matrix market format.\n");
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
	// "doc/man3/wls-check.3.gz"
	make_manpage(filename, &err);
	return err;
      }
  }
  
  bool verbose = cmd.get_logical("-v");
  WLS::integer_t nt;
  if (!cmd.get_integer("--nt",&nt))
    {
      if (verbose)
	{
	  std::cerr << "// WARNING, missing '--nt <numthreads>', default is 1." << std::endl;
	}
      nt = 1;
    }

  if (cmd.get_nargs() < 2)
    {
      std::cerr << "ERROR missing file" << std::endl;
      return 1;
    }

  const char * ifilename = cmd.get_arg(1);
  using real_t = double;
 
  WLS::integer_t * 	mat_begin = nullptr;
  WLS::integer_t * 	mat_idx = nullptr;
  real_t * 		mat_values = nullptr;
 

  WLS::integer_t 	mat_nrows;
  WLS::integer_t  	mat_ncols;
  WLS::integer_t  	mat_nnz;
  WLS::integer_t nmax_per_row, nmin_per_row;

  WLS::Input::matrix_market_t matrix_market;
  WLS::status_t status;

  status = WLS::Input::matrix_market_t::create(&matrix_market,ifilename);
  if (WLS::status_t::success != status)
    {
      fprintf(stderr,"status(=" iformat ")", (WLS::integer_t)status);
      return status;
    }
  
  status = matrix_market.dimension(&mat_nrows,
						 &mat_ncols,
						 &mat_nnz);
  if (WLS::status_t::success != status)
    {
      fprintf(stderr,"status(=" iformat ")", (WLS::integer_t)status);
      return status;
    }


  //
  // Allocate.
  //
  mat_begin 	= (WLS::integer_t*)calloc((mat_nrows+1),sizeof(WLS::integer_t));
  if (nullptr == mat_begin)
    {
      status = status_t::error_memory;
      goto state_error;
    }
   
  mat_idx	= (WLS::integer_t*)malloc(sizeof(WLS::integer_t)*(mat_nnz));
  if (nullptr == mat_idx)
    {
      status = status_t::error_memory;
      goto state_error;
    }
    
  mat_values 	= (real_t*)malloc(sizeof(real_t)*(mat_nnz));
  if (nullptr == mat_values)
    {
      status = status_t::error_memory;
      goto state_error;
    }



  status = matrix_market.nnz(WLS::mat_storage_t::row,
			     mat_nrows,
			     mat_ncols,
			     mat_nnz,
			     mat_begin);
  if (status)
    {
      fprintf(stderr,"error(=" iformat ")", (WLS::integer_t)status);
      return status;
    }
    
  WLS::partial_sums(mat_nrows,
		    mat_begin);

  status = matrix_market.create_compressed_sparse(WLS::mat_indexing_t::C,
						  WLS::mat_storage_t::row,
						  mat_nrows,
						  mat_ncols,
						  mat_nnz,
						  mat_begin,
						  mat_idx,
						  mat_values);
  if (status)
    {
      std::cerr << "matrix_market_read_cs failed, error code " << status << std::endl;
      return status;
    }

  nmax_per_row = mat_begin[1] - mat_begin[0];
  nmin_per_row = mat_begin[1] - mat_begin[0];
  for (WLS::integer_t irow = 1;irow < mat_nrows;++irow)
    {
      WLS::integer_t n_per_row = mat_begin[irow+1] - mat_begin[irow];
      nmax_per_row = std::max(nmax_per_row,n_per_row);	
      nmin_per_row = std::min(nmin_per_row,n_per_row);	
    }
   
  std::cout << "number of rows         : "
	    << mat_nrows
	    << std::endl;
  std::cout << "number of columns      : "
	    << mat_ncols
	    << std::endl;
  std::cout << "number of coefficients : "
	    << mat_nnz
	    << std::endl;
  std::cout << "nmin_per_row           : "
	    << nmin_per_row
	    << std::endl;    
  std::cout << "nmax_per_row           : "
	    << nmax_per_row
	    << std::endl;    
   
  free(mat_values);
  free(mat_idx);
  free(mat_begin);
  return 0;
  
 state_error:  
  {
    if (nullptr != mat_values)
      {
	free(mat_values);
	mat_values = nullptr;
      }

    if (nullptr != mat_idx)
      {
	free(mat_idx);
	mat_idx = nullptr;
      }

    if (nullptr != mat_begin)
      {
	free(mat_begin);
	mat_begin = nullptr;
      }
    return status;
  }
 
}
