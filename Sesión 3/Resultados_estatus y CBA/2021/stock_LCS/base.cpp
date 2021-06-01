  #include  <admodel.h>
  ofstream mcmc_report("mcmc.txt");
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <base.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nanos.allocate("nanos");
  nedades.allocate("nedades");
  edad_ini.allocate("edad_ini");
  delta_edad.allocate("delta_edad");
  ntallas.allocate("ntallas");
  indices.allocate(1,nanos,1,7,"indices");
  Tallas.allocate(1,ntallas,"Tallas");
  Cl.allocate(1,nanos,1,ntallas,"Cl");
  Nlcruceros.allocate(1,nanos,1,ntallas,"Nlcruceros");
  Wmed.allocate(1,nanos,1,ntallas,"Wmed");
  msex.allocate(1,ntallas,"msex");
  cvar.allocate(1,8,"cvar");
  nmus.allocate(1,2,"nmus");
  opt_qCru.allocate("opt_qCru");
  opt_devq.allocate(1,nanos,"opt_devq");
  opt_qCPUE.allocate("opt_qCPUE");
  opt_devqCPUE.allocate(1,nanos,"opt_devqCPUE");
  opt_Sel1.allocate("opt_Sel1");
  opt_Sel2.allocate(1,nanos,"opt_Sel2");
  opt_Sel3.allocate("opt_Sel3");
  opt_Sel4.allocate("opt_Sel4");
  opt_Sel5.allocate("opt_Sel5");
  parbiol.allocate(1,4,"parbiol");
  opt_VB1.allocate("opt_VB1");
  opt_VB2.allocate("opt_VB2");
  opt_VB3.allocate("opt_VB3");
  opt_VB4.allocate("opt_VB4");
  opt_VB5.allocate("opt_VB5");
  opt_Rmed.allocate("opt_Rmed");
  opt_devR.allocate("opt_devR");
  opt_devNo.allocate("opt_devNo");
  opt_F.allocate("opt_F");
  opt_M.allocate("opt_M");
  nanos_proy.allocate("nanos_proy");
  npbr.allocate("npbr");
  pR.allocate("pR");
  Fpbr.allocate(1,npbr,"Fpbr");
}

void model_parameters::initializationfunction(void)
{
  log_Rmed.set_initial_value(6);
  log_F.set_initial_value(-1.0);
  log_A50f_one.set_initial_value(1.09);
  log_Df_one.set_initial_value(0);
  log_A50cru.set_initial_value(1.09);
  log_Dcru.set_initial_value(0);
  log_sf.set_initial_value(0);
  log_sda.set_initial_value(0.40);
  log_sda2.set_initial_value(0.40);
  log_qCru.set_initial_value(0);
  dev_log_A50f.set_initial_value(0);
  dev_log_Df.set_initial_value(0);
  log_Lo.set_initial_value(2.9);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  log_A50f_one.allocate(0.7,1.6,opt_Sel1,"log_A50f_one");
  log_Df_one.allocate(0,0.7,opt_Sel1,"log_Df_one");
  dev_log_A50f.allocate(1,nanos,opt_Sel2,"dev_log_A50f");
  dev_log_Df.allocate(1,nanos,opt_Sel2,"dev_log_Df");
  log_muf.allocate(-9,1.24,opt_Sel3,"log_muf");
  log_sf.allocate(1,2,-10,3,opt_Sel3,"log_sf");
  log_A50cru.allocate(0.6,1.7,opt_Sel4,"log_A50cru");
  log_Dcru.allocate(0,0.7,opt_Sel4,"log_Dcru");
  log_Rmed.allocate(opt_Rmed,"log_Rmed");
  log_desv_No.allocate(1,nedades,-10,10,opt_devNo,"log_desv_No");
  log_desv_Rt.allocate(1,nanos-1,-10,10,opt_devR,"log_desv_Rt");
  log_F.allocate(1,nanos,-20,1.5,opt_F,"log_F");
  log_M.allocate(-3,1.5,opt_M,"log_M");
  log_qCru.allocate(opt_qCru,"log_qCru");
  log_qCPUE.allocate(opt_qCPUE,"log_qCPUE");
  devq.allocate(1,nanos-1,opt_devq,"devq");
  devqCPUE.allocate(1,nanos-1,opt_devqCPUE,"devqCPUE");
  log_Lo.allocate(2.8,3.2,opt_VB1,"log_Lo");
  log_cva.allocate(-3.9,-0.7,opt_VB2,"log_cva");
  log_sda.allocate(-2.9,1.2,opt_VB3,"log_sda");
  log_sda2.allocate(1,nedades,-2.9,1.2,opt_VB4,"log_sda2");
  log_k.allocate(-1.5,-0.5,opt_VB5,"log_k");
  ano.allocate(1,nanos,"ano");
  #ifndef NO_AD_INITIALIZE
    ano.initialize();
  #endif
  Desemb.allocate(1,nanos,"Desemb");
  #ifndef NO_AD_INITIALIZE
    Desemb.initialize();
  #endif
  Bcrucero.allocate(1,nanos,"Bcrucero");
  #ifndef NO_AD_INITIALIZE
    Bcrucero.initialize();
  #endif
  CPUE.allocate(1,nanos,"CPUE");
  #ifndef NO_AD_INITIALIZE
    CPUE.initialize();
  #endif
  cv1.allocate(1,nanos,"cv1");
  #ifndef NO_AD_INITIALIZE
    cv1.initialize();
  #endif
  cv2.allocate(1,nanos,"cv2");
  #ifndef NO_AD_INITIALIZE
    cv2.initialize();
  #endif
  cv3.allocate(1,nanos,"cv3");
  #ifndef NO_AD_INITIALIZE
    cv3.initialize();
  #endif
  Unos_edad.allocate(1,nedades,"Unos_edad");
  #ifndef NO_AD_INITIALIZE
    Unos_edad.initialize();
  #endif
  Unos_tallas.allocate(1,ntallas,"Unos_tallas");
  #ifndef NO_AD_INITIALIZE
    Unos_tallas.initialize();
  #endif
  Unos_ano.allocate(1,nanos,"Unos_ano");
  #ifndef NO_AD_INITIALIZE
    Unos_ano.initialize();
  #endif
  mu_edad.allocate(1,nedades,"mu_edad");
  #ifndef NO_AD_INITIALIZE
    mu_edad.initialize();
  #endif
  sigma_edad.allocate(1,nedades,"sigma_edad");
  #ifndef NO_AD_INITIALIZE
    sigma_edad.initialize();
  #endif
  Bcru.allocate(1,nanos,"Bcru");
  #ifndef NO_AD_INITIALIZE
    Bcru.initialize();
  #endif
  prior.allocate(1,10,"prior");
  #ifndef NO_AD_INITIALIZE
    prior.initialize();
  #endif
  Neq.allocate(1,nedades,"Neq");
  #ifndef NO_AD_INITIALIZE
    Neq.initialize();
  #endif
  Neqv.allocate(1,nedades,"Neqv");
  #ifndef NO_AD_INITIALIZE
    Neqv.initialize();
  #endif
  likeval.allocate(1,10,"likeval");
  #ifndef NO_AD_INITIALIZE
    likeval.initialize();
  #endif
  SDo.allocate(1,nanos,"SDo");
  #ifndef NO_AD_INITIALIZE
    SDo.initialize();
  #endif
  edades.allocate(1,nedades,"edades");
  #ifndef NO_AD_INITIALIZE
    edades.initialize();
  #endif
  Scru_1.allocate(1,nedades,"Scru_1");
  #ifndef NO_AD_INITIALIZE
    Scru_1.initialize();
  #endif
  Scru_2.allocate(1,nedades,"Scru_2");
  #ifndef NO_AD_INITIALIZE
    Scru_2.initialize();
  #endif
  log_A50f.allocate(1,nanos,"log_A50f");
  #ifndef NO_AD_INITIALIZE
    log_A50f.initialize();
  #endif
  log_Df.allocate(1,nanos,"log_Df");
  #ifndef NO_AD_INITIALIZE
    log_Df.initialize();
  #endif
  log_A50R2.allocate(1,nanos,"log_A50R2");
  #ifndef NO_AD_INITIALIZE
    log_A50R2.initialize();
  #endif
  log_DR2.allocate(1,nanos,"log_DR2");
  #ifndef NO_AD_INITIALIZE
    log_DR2.initialize();
  #endif
  qCru.allocate(1,nanos,"qCru");
  #ifndef NO_AD_INITIALIZE
    qCru.initialize();
  #endif
  qCPUE.allocate(1,nanos,"qCPUE");
  #ifndef NO_AD_INITIALIZE
    qCPUE.initialize();
  #endif
  Sflo.allocate(1,nanos,1,nedades,"Sflo");
  #ifndef NO_AD_INITIALIZE
    Sflo.initialize();
  #endif
  Scru.allocate(1,nanos,1,nedades,"Scru");
  #ifndef NO_AD_INITIALIZE
    Scru.initialize();
  #endif
  F.allocate(1,nanos,1,nedades,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Z.allocate(1,nanos,1,nedades,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  S.allocate(1,nanos,1,nedades,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  N.allocate(1,nanos,1,nedades,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  NM.allocate(1,nanos,1,nedades,"NM");
  #ifndef NO_AD_INITIALIZE
    NM.initialize();
  #endif
  Nv.allocate(1,nanos,1,nedades,"Nv");
  #ifndef NO_AD_INITIALIZE
    Nv.initialize();
  #endif
  Cedad.allocate(1,nanos,1,nedades,"Cedad");
  #ifndef NO_AD_INITIALIZE
    Cedad.initialize();
  #endif
  Prob_talla.allocate(1,nedades,1,ntallas,"Prob_talla");
  #ifndef NO_AD_INITIALIZE
    Prob_talla.initialize();
  #endif
  P1.allocate(1,nedades,1,ntallas,"P1");
  #ifndef NO_AD_INITIALIZE
    P1.initialize();
  #endif
  P2.allocate(1,nedades,1,ntallas,"P2");
  #ifndef NO_AD_INITIALIZE
    P2.initialize();
  #endif
  P3.allocate(1,nedades,1,ntallas,"P3");
  #ifndef NO_AD_INITIALIZE
    P3.initialize();
  #endif
  Cl_pred.allocate(1,nanos,1,ntallas,"Cl_pred");
  #ifndef NO_AD_INITIALIZE
    Cl_pred.initialize();
  #endif
  Nlcruceros_pred.allocate(1,nanos,1,ntallas,"Nlcruceros_pred");
  #ifndef NO_AD_INITIALIZE
    Nlcruceros_pred.initialize();
  #endif
  pobs.allocate(1,nanos,1,ntallas,"pobs");
  #ifndef NO_AD_INITIALIZE
    pobs.initialize();
  #endif
  ppred.allocate(1,nanos,1,ntallas,"ppred");
  #ifndef NO_AD_INITIALIZE
    ppred.initialize();
  #endif
  pobs_cru.allocate(1,nanos,1,ntallas,"pobs_cru");
  #ifndef NO_AD_INITIALIZE
    pobs_cru.initialize();
  #endif
  ppred_cru.allocate(1,nanos,1,ntallas,"ppred_cru");
  #ifndef NO_AD_INITIALIZE
    ppred_cru.initialize();
  #endif
  suma1.allocate("suma1");
  #ifndef NO_AD_INITIALIZE
  suma1.initialize();
  #endif
  suma2.allocate("suma2");
  #ifndef NO_AD_INITIALIZE
  suma2.initialize();
  #endif
  suma3.allocate("suma3");
  #ifndef NO_AD_INITIALIZE
  suma3.initialize();
  #endif
  suma4.allocate("suma4");
  #ifndef NO_AD_INITIALIZE
  suma4.initialize();
  #endif
  pStotf.allocate("pStotf");
  #ifndef NO_AD_INITIALIZE
  pStotf.initialize();
  #endif
  pSf.allocate("pSf");
  #ifndef NO_AD_INITIALIZE
  pSf.initialize();
  #endif
  penalty.allocate("penalty");
  #ifndef NO_AD_INITIALIZE
  penalty.initialize();
  #endif
  Linf.allocate("Linf");
  #ifndef NO_AD_INITIALIZE
  Linf.initialize();
  #endif
  k.allocate("k");
  #ifndef NO_AD_INITIALIZE
  k.initialize();
  #endif
  cv_edad.allocate("cv_edad");
  #ifndef NO_AD_INITIALIZE
  cv_edad.initialize();
  #endif
  sd_edad.allocate("sd_edad");
  #ifndef NO_AD_INITIALIZE
  sd_edad.initialize();
  #endif
  Lo.allocate("Lo");
  #ifndef NO_AD_INITIALIZE
  Lo.initialize();
  #endif
  M.allocate("M");
  #ifndef NO_AD_INITIALIZE
  M.initialize();
  #endif
  Nvplus.allocate("Nvplus");
  #ifndef NO_AD_INITIALIZE
  Nvplus.initialize();
  #endif
  Npplus.allocate("Npplus");
  #ifndef NO_AD_INITIALIZE
  Npplus.initialize();
  #endif
  Yp.allocate("Yp");
  #ifndef NO_AD_INITIALIZE
  Yp.initialize();
  #endif
  factor.allocate("factor");
  #ifndef NO_AD_INITIALIZE
  factor.initialize();
  #endif
  temp0.allocate(1,nedades,"temp0");
  #ifndef NO_AD_INITIALIZE
    temp0.initialize();
  #endif
  temp1.allocate(1,nedades,"temp1");
  #ifndef NO_AD_INITIALIZE
    temp1.initialize();
  #endif
  Fx.allocate("Fx");
  #ifndef NO_AD_INITIALIZE
  Fx.initialize();
  #endif
  Wedad.allocate(1,nedades,"Wedad");
  #ifndef NO_AD_INITIALIZE
    Wedad.initialize();
  #endif
  Np.allocate(1,nedades,"Np");
  #ifndef NO_AD_INITIALIZE
    Np.initialize();
  #endif
  NMp.allocate(1,nedades,"NMp");
  #ifndef NO_AD_INITIALIZE
    NMp.initialize();
  #endif
  Sp.allocate(1,nedades,"Sp");
  #ifndef NO_AD_INITIALIZE
    Sp.initialize();
  #endif
  Fp.allocate(1,nedades,"Fp");
  #ifndef NO_AD_INITIALIZE
    Fp.initialize();
  #endif
  Zp.allocate(1,nedades,"Zp");
  #ifndef NO_AD_INITIALIZE
    Zp.initialize();
  #endif
  Cap.allocate(1,ntallas,"Cap");
  #ifndef NO_AD_INITIALIZE
    Cap.initialize();
  #endif
  YTP.allocate(1,nanos_proy,1,npbr,"YTP");
  #ifndef NO_AD_INITIALIZE
    YTP.initialize();
  #endif
  Bp.allocate(1,nanos_proy,1,npbr,"Bp");
  #ifndef NO_AD_INITIALIZE
    Bp.initialize();
  #endif
  SDp.allocate(1,nanos_proy,1,npbr,"SDp");
  #ifndef NO_AD_INITIALIZE
    SDp.initialize();
  #endif
  Nvp.allocate(1,nedades,"Nvp");
  #ifndef NO_AD_INITIALIZE
    Nvp.initialize();
  #endif
  Sfp.allocate(1,nedades,"Sfp");
  #ifndef NO_AD_INITIALIZE
    Sfp.initialize();
  #endif
  SDvp.allocate(1,nanos_proy,"SDvp");
  #ifndef NO_AD_INITIALIZE
    SDvp.initialize();
  #endif
  CPUE_pred.allocate(1,nanos,"CPUE_pred");
  Bcrucero_pred.allocate(1,nanos,"Bcrucero_pred");
  Desemb_pred.allocate(1,nanos,"Desemb_pred");
  SD.allocate(1,nanos,"SD");
  RPR.allocate(1,nanos,"RPR");
  RPR2.allocate(1,nanos,"RPR2");
  Reclutas.allocate(1,nanos-1,"Reclutas");
  Flast.allocate("Flast");
  CTP.allocate(1,npbr,"CTP");
  RPRp.allocate(1,npbr,"RPRp");
  Btot.allocate(1,nanos,"Btot");
  Bv.allocate(1,nanos,"Bv");
  SDo2.allocate("SDo2");
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  ano=column(indices,1);// asigno la 1 columna de indices a "años"
  Bcrucero=column(indices,2);
  cv1=column(indices,3);
  CPUE=column(indices,4);
  cv2=column(indices,5);
  Desemb=column(indices,6);
  cv3=column(indices,7);
  Unos_edad=1;// lo uso en  operaciones matriciales con la edad
  Unos_ano=1;// lo uso en operaciones matriciales con el año
  Unos_tallas=1;// lo uso en operaciones matriciales con el año
  reporte_mcmc=0;
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{200,1000,3000,5000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-3,1e-5,1e-5,1e-6}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::userfunction(void)
{
  f =0.0;
  Eval_selectividad();
  Eval_mortalidades();
  Eval_abundancia();
  Eval_prob_talla_edad();
  Eval_capturas_predichas();
  Eval_deinteres();
  Eval_logverosim();
  Eval_funcion_objetivo();
  if (last_phase()){
  Eval_CTP();
  Eval_mcmc();}
}

void model_parameters::Eval_selectividad(void)
{
  int i;
  edades.fill_seqadd(edad_ini,delta_edad);
  log_A50f(1)=log_A50f_one;
  log_Df(1)=log_Df_one;
  for (int i=2;i<=nanos;i++){
  log_A50f(i)=log_A50f(i-1)+dev_log_A50f(i-1);
  log_Df(i)=log_Df(i-1)+dev_log_Df(i-1);}
  Scru_1=elem_div(Unos_edad,(1+exp(-1.0*log(19)*(edades-exp(log_A50cru))/exp(log_Dcru))));
  for (i=1;i<=nanos;i++)
  {Sflo(i)=elem_div(Unos_edad,(1+mfexp(-1.0*log(19)*(edades-mfexp(log_A50f(i)))/mfexp(log_Df(i)))));
  Scru(i)=Scru_1;
  }
  if(opt_Sel3>0){// selectividad doble_normal unica 
    for (i=1;i<=nanos;i++){
      Sflo(i)=mfexp(-1/(2*square(exp(log_sf(1))))*square((edades-exp(log_muf))));
      for (int j=1;j<=nedades;j++){
       if(edades(j)>exp(log_muf)){
       Sflo(i,j)=mfexp(-1/(2*square(exp(log_sf(2))))*square(-1.*(edades(j)-exp(log_muf))));}
    }}}
  if(opt_Sel5>0)// selectividad 1.0
      {Scru=1;}
}

void model_parameters::Eval_mortalidades(void)
{
  M=parbiol(4);
  if (active(log_M)){M=mfexp(log_M);}
  F=elem_prod(outer_prod(mfexp(log_F),Unos_edad),Sflo);
  Z=F+M;
  S=mfexp(-1.0*Z);
  Flast=mfexp(log_F(nanos));
}

void model_parameters::Eval_abundancia(void)
{
 int i, j;
  Neq(1)=mfexp(log_Rmed);
  Neqv(1)=mfexp(log_Rmed);
  for (j=2;j<=nedades;j++)
  { Neq(j)=Neq(j-1)*exp(-1.*Z(1,j-1));// stock eq
    Neqv(j)=Neqv(j-1)*exp(-1.*M);} // stock virginal en eq
    Neq(nedades)=Neq(nedades)/(1-exp(-1.*Z(1,nedades)));// grupo plus
    Neqv(nedades)=Neqv(nedades)/(1-exp(-1.*M));// grupo plus
   N(1)=mfexp(log(Neq)+log_desv_No);// Composición de ededades poblacion inicial
   Reclutas=mfexp(log_Rmed+log_desv_Rt);
  for (i=2;i<=nanos;i++)
  {N(i,1)=Reclutas(i-1);}
  for (i=1;i<nanos;i++)
  {N(i+1)(2,nedades)=++elem_prod(N(i)(1,nedades-1),S(i)(1,nedades-1));
   N(i+1,nedades)+=N(i,nedades)*S(i,nedades);
  }
}

void model_parameters::Eval_prob_talla_edad(void)
{
  Linf=parbiol(1);
  k=parbiol(2);
  Lo=parbiol(3);
  if (active(log_Lo)){Lo=mfexp(log_Lo);}
  if (active(log_k)){k=mfexp(log_k);}
 int i, j;
  mu_edad(1)=Lo;
  for (i=2;i<=nedades;i++){
   mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);
   }
    cv_edad=mfexp(log_cva);
    sigma_edad=cv_edad*mu_edad;  // proporcional con la talla por defecto
  if (active(log_sda)){
    sigma_edad=mfexp(log_sda);}// constante
  if (active(log_sda2)){// independiente
   sigma_edad=mfexp(log_sda2);
    }
  for (i=1;i<=nedades;i++){
    P1(i)=(Tallas-mu_edad(i))/sigma_edad(i);
  for (j=1;j<=ntallas;j++){
    P2(i,j)=cumd_norm(P1(i,j));}}
  for (i=1;i<=nedades;i++){
     for (j=2;j<=ntallas;j++){
       P3(i,j)=P2(i,j)-P2(i,j-1);}}
  Prob_talla=elem_div(P3+1e-16,outer_prod(rowsum(P3+1e-16),Unos_tallas));
}

void model_parameters::Eval_capturas_predichas(void)
{
  Cedad=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));
  Cl_pred=Cedad*Prob_talla;
  NM=elem_div(elem_prod(N,1-S),Z);   
  Nlcruceros_pred=elem_prod(NM,Scru)*Prob_talla;
  pobs=elem_div(Cl,outer_prod(rowsum(Cl),Unos_tallas)+1e-5);
  ppred=elem_div(Cl_pred,outer_prod(rowsum(Cl_pred),Unos_tallas));
  pobs_cru=elem_div(Nlcruceros,outer_prod(rowsum(Nlcruceros),Unos_tallas)+1e-5);
  ppred_cru=elem_div(Nlcruceros_pred,outer_prod(rowsum(Nlcruceros_pred),Unos_tallas)+1e-5);
  Desemb_pred=rowsum((elem_prod(Cl_pred,Wmed)));
}

void model_parameters::Eval_deinteres(void)
{
  Nv=N;// solo para empezar los calculos
 for (int i=1;i<nanos;i++)
  {
      Nv(i+1)(2,nedades)=++Nv(i)(1,nedades-1)*exp(-1.0*M);
      Nv(i+1,nedades)+=Nv(i,nedades)*exp(-1.0*M);}
  Btot=rowsum((elem_prod(N*Prob_talla,Wmed)));
  Bcru=rowsum((elem_prod(Nlcruceros_pred,Wmed)));;// biomasas al crucero 
  Bv=rowsum(elem_prod(elem_prod(NM,Sflo)*Prob_talla,Wmed));
   qCru(1)=mfexp(log_qCru);
   qCPUE(1)=mfexp(log_qCPUE);
   for (int i=2;i<=nanos;i++)
   {qCru(i)=qCru(i-1)*mfexp(devq(i-1));
    qCPUE(i)=qCPUE(i-1)*mfexp(devqCPUE(i-1));}
   Bcrucero_pred=elem_prod(qCru,Bcru);
  SD=rowsum(elem_prod(elem_prod(elem_prod(N,exp(-0.67*Z))*Prob_talla,Wmed),outer_prod(Unos_ano,msex)));// desovantes al 1 agosto
  SDo=rowsum(elem_prod(elem_prod(Nv*exp(-0.67*M)*Prob_talla,Wmed),outer_prod(Unos_ano,msex)));// sin pesca al 1 agosto
  SDo2=sum(elem_prod(elem_prod((Neqv*exp(-0.67*M))*Prob_talla,msex),Wmed(nanos)));// virginal al 1 agosto
  RPR=elem_div(SD,SDo);
  RPR2=SD/SDo2;
  CPUE_pred=elem_prod(qCPUE,Bv);
}

void model_parameters::Eval_logverosim(void)
{
  int i;
  suma1=0; suma2=0; 
  for (i=1;i<=nanos;i++)
  {
   if (Bcrucero(i)>0){
    suma1+=square((log(Bcrucero(i))-log(Bcrucero_pred(i)))/cv1(i));}
   if (CPUE(i)>0){
    suma2+=square((log(CPUE(i))-log(CPUE_pred(i)))/cv2(i));}
  }
}

void model_parameters::Eval_funcion_objetivo(void)
{
  int i,j;
  likeval(1)=0.5*suma1;//Cruceros
  likeval(2)=0.5*suma2;//CPUE
  likeval(3)=0.5*norm2(elem_div(log(Desemb)-log(Desemb_pred),cv3));//Desembarques  
  likeval(4)=-1.*nmus(1)*sum(elem_prod(pobs,log(ppred)));
  likeval(5)=-1.*nmus(2)*sum(elem_prod(pobs_cru,log(ppred_cru)));
  prior(1)=0.5*norm2((log(row(N,1))-log(Neq))/cvar(1));
  prior(2)=0.5*norm2(log_desv_Rt/cvar(2));
  if (active(log_k)){// si estima k
  prior(3)=0.5*square((log_k-log(parbiol(2)))/cvar(3));}
  if (active(log_qCru)){
  prior(4)=0.5*square(log_qCru/cvar(4));}
  if(max(opt_Sel2)>0){
  prior(5)=0.5*norm2(dev_log_A50f/cvar(5));}
  if (active(log_M)){
  prior(6)=0.5*square((log_M-parbiol(4))/cvar(6));}
  if (max(opt_devq)>0){
  prior(7)=0.5*norm2(devq/cvar(7));}
  if (max(opt_devqCPUE)>0){
  prior(8)=0.5*norm2(devqCPUE/cvar(8));}
   f=sum(likeval)+sum(prior);
}

void model_parameters::Eval_CTP(void)
{
 Sfp=Sflo(nanos);//elem_div(Unos_edad,(1+mfexp(-1.0*log(19)*(edades-mfexp(log_A50f_one))/mfexp(log_Df_one))));
 for (int j=1;j<=npbr;j++){ // pbr
  Np=N(nanos);
  NMp=NM(nanos);
  Sp=S(nanos);
 for (int i=1;i<=nanos_proy;i++){
  Npplus=Np(nedades)*Sp(nedades);
  Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
  Np(nedades)+=Npplus;
  Np(1)= pR*mfexp(log_Rmed);//mean(Reclutas); // proyecta con R promedi
  Bp(i,j)=sum(elem_prod(Np*Prob_talla,Wmed(nanos)));
   Fp = Sfp*Fpbr(j);
   Zp = Fp + M;
   Sp = mfexp(-1.0*(Zp));
  NMp=elem_prod(Np,pow(Sp,0.67));
 // SDp(i,j)=sum(elem_prod(elem_prod(Np,exp(-0.67*Zp))*Prob_talla,elem_prod(Wmed(nanos),msex)));
  SDp(i,j)=sum(elem_prod(NMp*Prob_talla,elem_prod(Wmed(nanos),msex)));
  Cap=elem_prod(elem_div(Fp,Zp),elem_prod(Np,(1-Sp)))*Prob_talla;
  YTP(i,j)=sum(elem_prod(Cap,Wmed(nanos)));
  }
  CTP(j)=YTP(1,j);
  }
 // Rutina para la estimación de RPR
 Nvp=Nv(nanos);// toma la ultima estimación
 for (int i=1;i<=nanos_proy;i++)
  {
      Nvplus=Nvp(nedades)*exp(-1.0*M);
      Nvp(2,nedades)=++Nvp(1,nedades-1)*exp(-1.0*M);
      Nvp(nedades)+=Nvplus;
      Nvp(1)= mfexp(log_Rmed);//mean(Reclutas);
      SDvp(i)=sum(elem_prod((Nvp*exp(-0.67*M))*Prob_talla,elem_prod(Wmed(nanos),msex)));
  }
 // mido el riesgo de RPR solo para el año proyectado de interés
 for (int i=1;i<=npbr;i++)
  {
  RPRp(i)=SDp(1,i)/SDvp(i);//
  }
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  report << "Bcru_obs" << endl;
  report << Bcrucero << endl;
  report << "Bcru_pred" << endl;
  report << Bcrucero_pred << endl;
  report << "Y_obs" << endl;
  report << Desemb << endl;
  report << "Y_pred" << endl;
  report << Desemb_pred << endl;
  report << "cpue_obs" << endl;
  report << CPUE << endl;
  report << "cpue_pred" << endl;
  report << CPUE_pred << endl;
  report << "Btot" << endl;
  report << Btot << endl;
  report << "BD" << endl;
  report << SD << endl;
  report << "BV" << endl;
  report << Bv << endl;
  report << "RPR" << endl;
  report << RPR2 << endl;
  report << "N " << endl;
  report << N << endl;
  report << "Neq " << endl;
  report << Neq << endl;
  report << "F" << endl;
  report << F << endl;
  report << "Sflo" << endl;
  report << Sflo << endl;
  report << "pf_obs " << endl;
  report << pobs << endl;
  report << "pf_pred " << endl;
  report << ppred << endl;
  report << "pcru_obs " << endl;
  report << pobs_cru << endl;
  report << "pcru_pred " << endl;
  report << ppred_cru << endl;
  report << "Scru" << endl;
  report << Scru << endl;
  report << "Lf_obs_pred" << endl;
  report << Tallas*trans(pobs)<< endl;
  report << Tallas*trans(ppred)<< endl;
  report << "Lc_obs_pred" << endl;
  report << Tallas*trans(pobs_cru)<< endl;
  report << Tallas*trans(ppred_cru)<< endl;
  report << "log_like " << endl;
  report << likeval << endl;
  report << "Crecimiento " << endl;
  report <<"Loo    k    Lo "<<endl;
  report << Linf <<" "<< k <<"  "<< Lo <<endl;
  report << "edades " << endl;
  report << edades << endl;
  report << "talla_edad " << endl;
  report << mu_edad << endl;
  report << "sigma_edad " << endl;
  report << sigma_edad << endl;
  report << "_____________________________ " << endl;
  report << "Capturabilidad " << endl;
  report <<"qcru"<<endl;
  report << qCru<<endl;
  report <<"BD0"<<endl;
  report << SDo2<<endl;
  report << "_____________________________ " << endl;
  report << "matriz_p" << endl;
  report << Prob_talla << endl;
  report << "_____________________________ " << endl;
  report << "BT_proy" << endl;
  report << Bp << endl;
  report << "Y_proy" << endl;
  report << YTP << endl;
  report << "BD_proy"<< endl;
  report << SDp<< endl;
  report << "_____________________________ " << endl;
  report << "M" << endl;
  report << M <<endl;
  report << "_____________________________ " << endl;
  report << "Yedad" << endl;
  report << Cedad <<endl;
  report << "RPR proyectado" << endl;
  report << RPRp << endl;
  report << SDp(nanos_proy)/SDo2 << endl;
  report << "FPBR" << endl;
  report << Fpbr << endl;
}

void model_parameters::Eval_mcmc(void)
{
  if(reporte_mcmc == 0)
  mcmc_report<<"Bcru CTP1 CTP2 CTP3 CTP4 CTP5 BDp1_fin BDp2_fin BDp3_fin BDp4_fin BDp5_fin"<<endl;
  mcmc_report<<Bcrucero_pred(nanos)<<" "<<YTP(1,1)<<" "<<YTP(1,2)<<" "<<YTP(1,3)<<" "<<YTP(1,4)<<" "<<YTP(1,5)<<
     " "<<SDp(nanos_proy,1)<<" "<<SDp(nanos_proy,2)<<" "<<SDp(nanos_proy,3)<<
     " "<<SDp(nanos_proy,4)<<" "<<SDp(nanos_proy,5)<<endl;
  reporte_mcmc++;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize=300000; // 
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(30000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(50000000);
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
