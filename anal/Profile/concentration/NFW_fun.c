struct data {
   size_t nbin;
   double * r;
   double * rho;
   double * sigma;
 };
 
 int NFW_f (const gsl_vector * par, void *data, 
		 gsl_vector * f)
 {
   size_t nbin = ((struct data *)data)->nbin;
   double *r=((struct data *)data)->r;
   double *rho = ((struct data *)data)->rho;
   double *sigma = ((struct data *) data)->sigma;
 
   double rhos = gsl_vector_get (par, 0);
   double rs = gsl_vector_get (par, 1);
 
   size_t i;
 
   for (i = 0; i < nbin; i++)
	 {
	   /* Model Yi = Rhos /rt/(1+rt)^2 */
	   double rt = r[i]/rs;
	   double Yi = rhos/rt/(1+rt)/(1+rt);
	   gsl_vector_set (f, i, (Yi - rho[i])/sigma[i]);
	 }
 
   return GSL_SUCCESS;
 }
 
 int NFW_df (const gsl_vector * par, void *data, 
		  gsl_matrix * J)
 {
   size_t nbin = ((struct data *)data)->nbin;
   double *r=((struct data *)data)->r;
   double *sigma = ((struct data *) data)->sigma;
 
   double rhos = gsl_vector_get (par, 0);
   double rs = gsl_vector_get (par, 1);
 
   size_t i;
 
   for (i = 0; i < nbin; i++)
	 {
	   /* Jacobian matrix J(i,j) = dfi / dxj, */
	   /* where fi = (Yi - yi)/sigma[i],      */
	   /*       Yi = Rhos /rt/(1+rt)^2  */
	   /* and the xj are the parameters (rhos,rs) */
	   double rt = r[i]/rs;
	   double s = sigma[i];
	   double den=1.0/rt/(1+rt)/(1+rt);
	   gsl_matrix_set (J, i, 0, den/s); 
	   gsl_matrix_set (J, i, 1, rhos/rs*(3-2./(1+rt))*den/s);
	 }
   return GSL_SUCCESS;
 }
 
 int NFW_fdf (const gsl_vector * x, void *data,
		   gsl_vector * f, gsl_matrix * J)
 {
   NFW_f (x, data, f);
   NFW_df (x, data, J);
 
   return GSL_SUCCESS;
 }

