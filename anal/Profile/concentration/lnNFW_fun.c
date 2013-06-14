struct data {
   size_t nbin;
   double * r;
   double * rho;//log(rho)
   double * sigma;//error for log(rho)
 };
 
 int lnNFW_f (const gsl_vector * par, void *data, 
		 gsl_vector * f)
 {
   size_t nbin = ((struct data *)data)->nbin;
   double *r=((struct data *)data)->r;
   double *Lrho = ((struct data *)data)->rho;//log(rho)
   double *Lsigma = ((struct data *) data)->sigma;//D(log(rho))
 
   double Lrhos = gsl_vector_get (par, 0);
   double rs = gsl_vector_get (par, 1);
 
   size_t i;
 
   for (i = 0; i < nbin; i++)
	 {
	   /* Model Yi = Rhos /rt/(1+rt)^2 */
	   double rt = r[i]/rs;
	   double Yi = Lrhos-log(rt)-2.*log(1+rt);
	   gsl_vector_set (f, i, (Yi - Lrho[i])/Lsigma[i]);
	 }
 
   return GSL_SUCCESS;
 }
 
 int lnNFW_df (const gsl_vector * par, void *data, 
		  gsl_matrix * J)
 {
   size_t nbin = ((struct data *)data)->nbin;
   double *r=((struct data *)data)->r;
   double *Lsigma = ((struct data *) data)->sigma;
 
   double rs = gsl_vector_get (par, 1);
 
   size_t i;
 
   for (i = 0; i < nbin; i++)
	 {
	   /* Jacobian matrix J(i,j) = dfi / dxj, */
	   /* where fi = (Yi - yi)/sigma[i],      */
	   /*   Yi=log(Rhoi),  Rhoi= Rhos /rt/(1+rt)^2  */
	   /* and the xj are the parameters (Lrhos,rs) */
	   double rt = r[i]/rs;
	   double s = Lsigma[i];
	   gsl_matrix_set (J, i, 0, 1./s); 
	   gsl_matrix_set (J, i, 1, (1.+2.*rt/(1+rt))/rs/s);
	 }
   return GSL_SUCCESS;
 }
 
 int lnNFW_fdf (const gsl_vector * x, void *data,
		   gsl_vector * f, gsl_matrix * J)
 {
   lnNFW_f (x, data, f);
   lnNFW_df (x, data, J);
 
   return GSL_SUCCESS;
 }

