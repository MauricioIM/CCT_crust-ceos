GLOBALS_SECTION
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
 ofstream mcmc_report("mcmc2.csv");

TOP_OF_MAIN_SECTION
  arrmblsize=300000; // 
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(30000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(50000000);
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);

DATA_SECTION
  init_int nanos  
  init_int nedades
  init_number edad_ini
  init_number delta_edad
  init_int ntallas
  init_matrix indices(1,nanos,1,7)
  init_vector Tallas(1,ntallas)
  init_matrix Cl(1,nanos,1,ntallas)
  init_matrix Nlcruceros(1,nanos,1,ntallas)
  init_matrix Wmed(1,nanos,1,ntallas)
  init_vector msex(1,ntallas)
  int reporte_mcmc

//!! ad_comm::change_datafile_name("lacnorte_2016.ctl");

  init_vector cvar(1,8)
  init_vector nmus(1,2)
  init_int    opt_qCru
  init_ivector opt_devq(1,nanos)
  init_int    opt_qCPUE
  init_ivector opt_devqCPUE(1,nanos)
  init_int    opt_Sel1
  init_ivector opt_Sel2(1,nanos)
  init_int    opt_Sel3
  init_int    opt_Sel4
  init_int    opt_Sel5

  init_vector parbiol(1,4)
  init_int    opt_VB1
  init_int    opt_VB2
  init_int    opt_VB3
  init_int    opt_VB4
  init_int    opt_VB5
  init_int    opt_Rmed
  init_int    opt_devR
  init_int    opt_devNo
  init_int    opt_F
  init_int    opt_M
  init_int    nanos_proy
  init_int    npbr
  init_number pR
  init_vector Fpbr(1,npbr)


INITIALIZATION_SECTION
// defino un valor inicial de log_reclutamiento promedio (factor de escala)
  log_Rmed         6
  log_F            -1.0
  log_A50f_one     1.09
  log_Df_one       0
  log_A50cru       1.09
  log_Dcru         0
  log_sf            0
  log_sda           0.40
  log_sda2          0.40
  log_qCru          0
  dev_log_A50f      0
  dev_log_Df        0
  log_Lo            2.9


PARAMETER_SECTION

// parametros selectividad
 init_bounded_number log_A50f_one(0.7,1.6,opt_Sel1)  
 init_bounded_number log_Df_one(0,0.7,opt_Sel1)

// desvios de los parametros selectividad flota
 init_number_vector dev_log_A50f(1,nanos,opt_Sel2)  
 init_number_vector dev_log_Df(1,nanos,opt_Sel2)

// parametros selectividad doble normal
 init_bounded_number log_muf(-9,1.24,opt_Sel3)
 init_bounded_vector log_sf(1,2,-10,3,opt_Sel3)

// parametros selectividad unica Cruceros
 init_bounded_number log_A50cru(0.6,1.7,opt_Sel4)  
 init_bounded_number log_Dcru(0,0.7,opt_Sel4)

// Recruits and F mortalities
// parametros reclutamientos y mortalidades)
 init_number log_Rmed(opt_Rmed)
 init_bounded_vector log_desv_No(1,nedades,-10,10,opt_devNo)
 init_bounded_vector log_desv_Rt(1,nanos-1,-10,10,opt_devR)
 init_bounded_vector log_F(1,nanos,-20,1.5,opt_F) // log  mortalidad por pesca por flota
 init_bounded_number log_M(-3,1.5,opt_M)

// capturabilidades
 init_number log_qCru(opt_qCru)
 init_number log_qCPUE(opt_qCPUE)
 init_number_vector devq(1,nanos-1,opt_devq)
 init_number_vector devqCPUE(1,nanos-1,opt_devqCPUE)

// crecimiento
 init_bounded_number log_Lo(2.8,3.2,opt_VB1)
 init_bounded_number log_cva(-3.9,-0.7,opt_VB2)
 init_bounded_number log_sda(-2.9,1.2,opt_VB3)
 init_bounded_vector log_sda2(1,nedades,-2.9,1.2,opt_VB4)
 init_bounded_number log_k(-1.5,-0.5,opt_VB5)

//Defino las variables de estado 
 vector ano(1,nanos)
 vector Desemb(1,nanos)
 //vector Desemb_pred(1,nanos);
 vector Bcrucero(1,nanos)
 //vector Bcrucero_pred(1,nanos);
 
 vector CPUE(1,nanos)
 vector cv1(1,nanos)
 vector cv2(1,nanos)
 vector cv3(1,nanos)
 vector Unos_edad(1,nedades)
 vector Unos_tallas(1,ntallas)
 vector Unos_ano(1,nanos)
 vector mu_edad(1,nedades);
 vector sigma_edad(1,nedades);
 vector Bcru(1,nanos);
 vector prior(1,10);
// vector CPUE_pred(1,nanos);
 vector Neq(1,nedades);
 vector Neqv(1,nedades);
 vector likeval(1,10);
 vector SDo(1,nanos);
 number SDo2;

 vector edades(1,nedades)
 vector Scru_1(1,nedades);
 vector Scru_2(1,nedades);
 vector log_A50f(1,nanos)
 vector log_Df(1,nanos)
 vector log_A50R2(1,nanos)
 vector log_DR2(1,nanos)
 vector qCru(1,nanos)
 vector qCPUE(1,nanos)

 matrix Sflo(1,nanos,1,nedades)
 matrix Scru(1,nanos,1,nedades)
 matrix F(1,nanos,1,nedades)
 matrix Z(1,nanos,1,nedades)
 matrix S(1,nanos,1,nedades)
 matrix N(1,nanos,1,nedades)
 matrix NM(1,nanos,1,nedades)
 matrix Nv(1,nanos,1,nedades)
 matrix Cedad(1,nanos,1,nedades)
 matrix Prob_talla(1,nedades,1,ntallas)
 matrix P1(1,nedades,1,ntallas)
 matrix P2(1,nedades,1,ntallas)
 matrix P3(1,nedades,1,ntallas)
 matrix Cl_pred(1,nanos,1,ntallas)
 matrix Nlcruceros_pred(1,nanos,1,ntallas)
 matrix pobs(1,nanos,1,ntallas)
 matrix ppred(1,nanos,1,ntallas)
 matrix pobs_cru(1,nanos,1,ntallas)
 matrix ppred_cru(1,nanos,1,ntallas)

 number suma1
 number suma2
 number suma3
 number suma4
 number pStotf
 number pSf
 number penalty

 number Linf
 number k
 number cv_edad
 number sd_edad
 number Lo
 number M
 number Nvplus
 number Npplus


// los arreglos usados en la proyeccion

  number Yp
  number factor
  vector temp0(1,nedades)
  vector temp1(1,nedades)
  number Fx
  vector Wedad(1,nedades)



 vector Np(1,nedades)
 vector NMp(1,nedades)
 vector Sp(1,nedades)
 vector Fp(1,nedades)
 vector Zp(1,nedades)
 vector Cap(1,ntallas)

 matrix YTP(1,nanos_proy,1,npbr)
 matrix Bp(1,nanos_proy,1,npbr)
 //matrix SDp(1,nanos_proy,1,npbr)

 vector Nvp(1,nedades)
 vector Sfp(1,nedades)
 vector SDvp(1,nanos_proy)

 sdreport_vector CPUE_pred(1,nanos) //
 sdreport_vector Bcrucero_pred(1,nanos) //
 sdreport_vector Desemb_pred(1,nanos) //
 sdreport_vector SD(1,nanos) // 
 sdreport_vector RPR(1,nanos) // 
 sdreport_vector RPR2(1,nanos) // 
 sdreport_vector Reclutas(1,nanos-1)
 sdreport_number Flast
 sdreport_vector CTP(1,npbr)
 sdreport_vector RPRp(1,npbr)
 sdreport_vector Btot(1,nanos)
 sdreport_vector Bv(1,nanos)
 sdreport_matrix SDp(1,nanos_proy,1,npbr)

  objective_function_value f


PRELIMINARY_CALCS_SECTION
// leo la matriz de indices

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

RUNTIME_SECTION
  maximum_function_evaluations 200,1000,3000,5000
  convergence_criteria  1e-3,1e-5,1e-5,1e-6


PROCEDURE_SECTION
// para comentar mas de una lina  /*.........*/

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



FUNCTION Eval_selectividad
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

// parametros selectividad doble normal

  if(opt_Sel3>0){// selectividad doble_normal unica 
    for (i=1;i<=nanos;i++){
      Sflo(i)=mfexp(-1/(2*square(exp(log_sf(1))))*square((edades-exp(log_muf))));
      for (int j=1;j<=nedades;j++){
       if(edades(j)>exp(log_muf)){
       Sflo(i,j)=mfexp(-1/(2*square(exp(log_sf(2))))*square(-1.*(edades(j)-exp(log_muf))));}
    }}}
     
  if(opt_Sel5>0)// selectividad 1.0
      {Scru=1;}


FUNCTION Eval_mortalidades


  M=parbiol(4);

  if (active(log_M)){M=mfexp(log_M);}
  F=elem_prod(outer_prod(mfexp(log_F),Unos_edad),Sflo);
  Z=F+M;
  S=mfexp(-1.0*Z);

  Flast=mfexp(log_F(nanos));

FUNCTION Eval_abundancia
 int i, j;

// genero una estructura inicial de equilibrio como referencia para el primer año
  Neq(1)=mfexp(log_Rmed);
  Neqv(1)=mfexp(log_Rmed);

  for (j=2;j<=nedades;j++)
  { Neq(j)=Neq(j-1)*exp(-1.*Z(1,j-1));// stock eq
    Neqv(j)=Neqv(j-1)*exp(-1.*M);} // stock virginal en eq

    Neq(nedades)=Neq(nedades)/(1-exp(-1.*Z(1,nedades)));// grupo plus
    Neqv(nedades)=Neqv(nedades)/(1-exp(-1.*M));// grupo plus

   N(1)=mfexp(log(Neq)+log_desv_No);// Composición de ededades poblacion inicial

   Reclutas=mfexp(log_Rmed+log_desv_Rt);

// luego considero los reclutas anuales a la edad 2
  for (i=2;i<=nanos;i++)
  {N(i,1)=Reclutas(i-1);}

// se estima la sobrevivencia por edad(a+1) y año(t+1)
  for (i=1;i<nanos;i++)
  {N(i+1)(2,nedades)=++elem_prod(N(i)(1,nedades-1),S(i)(1,nedades-1));
   N(i+1,nedades)+=N(i,nedades)*S(i,nedades);
  }


FUNCTION Eval_prob_talla_edad

  Linf=parbiol(1);
  k=parbiol(2);
  Lo=parbiol(3);

  if (active(log_Lo)){Lo=mfexp(log_Lo);}
  if (active(log_k)){k=mfexp(log_k);}

 int i, j;
 
// genero una clave edad-talla para otros calculos. Se modela desde L(1)
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

//---------------------------------------------------------------
  for (i=1;i<=nedades;i++){
    P1(i)=(Tallas-mu_edad(i))/sigma_edad(i);

  for (j=1;j<=ntallas;j++){
    P2(i,j)=cumd_norm(P1(i,j));}}
   
  for (i=1;i<=nedades;i++){
     for (j=2;j<=ntallas;j++){
       P3(i,j)=P2(i,j)-P2(i,j-1);}}

  Prob_talla=elem_div(P3+1e-16,outer_prod(rowsum(P3+1e-16),Unos_tallas));


FUNCTION Eval_capturas_predichas

// matrices de capturas predichas por edad y año
  Cedad=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));

// matrices de capturas predichas por talla y año
  Cl_pred=Cedad*Prob_talla;

// matrices de cruceros predichas por talla y año
  NM=elem_div(elem_prod(N,1-S),Z);   
  Nlcruceros_pred=elem_prod(NM,Scru)*Prob_talla;
 
// matrices de proporcion de capturas por talla y año
  pobs=elem_div(Cl,outer_prod(rowsum(Cl),Unos_tallas)+1e-5);
  ppred=elem_div(Cl_pred,outer_prod(rowsum(Cl_pred),Unos_tallas));

// Cruceros
  pobs_cru=elem_div(Nlcruceros,outer_prod(rowsum(Nlcruceros),Unos_tallas)+1e-5);
  ppred_cru=elem_div(Nlcruceros_pred,outer_prod(rowsum(Nlcruceros_pred),Unos_tallas)+1e-5);

// vectores de desembarques predichos por año
  Desemb_pred=rowsum((elem_prod(Cl_pred,Wmed)));


FUNCTION Eval_deinteres


// Rutina para calcular RPR
  Nv=N;// solo para empezar los calculos

 for (int i=1;i<nanos;i++)
  {
      Nv(i+1)(2,nedades)=++Nv(i)(1,nedades-1)*exp(-1.0*M);
      Nv(i+1,nedades)+=Nv(i,nedades)*exp(-1.0*M);}

  Btot=rowsum((elem_prod(N*Prob_talla,Wmed)));
  Bcru=rowsum((elem_prod(Nlcruceros_pred,Wmed)));;// biomasas al crucero 
  Bv=rowsum(elem_prod(elem_prod(NM,Sflo)*Prob_talla,Wmed));
  

// Estimacion de B.cruceros

   qCru(1)=mfexp(log_qCru);
   qCPUE(1)=mfexp(log_qCPUE);

   for (int i=2;i<=nanos;i++)
   {qCru(i)=qCru(i-1)*mfexp(devq(i-1));
    qCPUE(i)=qCPUE(i-1)*mfexp(devqCPUE(i-1));}

   Bcrucero_pred=elem_prod(qCru,Bcru);


  SD=rowsum(elem_prod(elem_prod(elem_prod(N,exp(-0.67*Z))*Prob_talla,Wmed),outer_prod(Unos_ano,msex)));// desovantes al 1 agosto
  SDo=rowsum(elem_prod(elem_prod(Nv*exp(-0.67*M)*Prob_talla,Wmed),outer_prod(Unos_ano,msex)));// sin pesca al 1 agosto

//  SDo2=sum(elem_prod(elem_prod((Neqv*exp(-0.67*M))*Prob_talla,msex),colsum(Wmed)/nanos));// virginal al 1 agosto
  SDo2=sum(elem_prod(elem_prod((Neqv*exp(-0.67*M))*Prob_talla,msex),Wmed(nanos)));// virginal al 1 agosto


  RPR=elem_div(SD,SDo);
  RPR2=SD/SDo2;

  CPUE_pred=elem_prod(qCPUE,Bv);


FUNCTION Eval_logverosim
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
  int i;

  suma1=0; suma2=0; 

  for (i=1;i<=nanos;i++)
  {
   if (Bcrucero(i)>0){
    suma1+=square((log(Bcrucero(i))-log(Bcrucero_pred(i)))/cv1(i));}
   if (CPUE(i)>0){
    suma2+=square((log(CPUE(i))-log(CPUE_pred(i)))/cv2(i));}
  }


FUNCTION Eval_funcion_objetivo

  int i,j;

// se calcula la F.O. como la suma de las -logver
// lognormalgraf
  likeval(1)=0.5*suma1;//Cruceros
  likeval(2)=0.5*suma2;//CPUE
  likeval(3)=0.5*norm2(elem_div(log(Desemb)-log(Desemb_pred),cv3));//Desembarques  

// multinomial flota
  likeval(4)=-1.*nmus(1)*sum(elem_prod(pobs,log(ppred)));

// multinomial cruceros
  likeval(5)=-1.*nmus(2)*sum(elem_prod(pobs_cru,log(ppred_cru)));


// Priors
// lognormal Ninicial y Reclutas
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


FUNCTION Eval_CTP

 Sfp=Sflo(nanos);//elem_div(Unos_edad,(1+mfexp(-1.0*log(19)*(edades-mfexp(log_A50f_one))/mfexp(log_Df_one))));

 for (int j=1;j<=npbr;j++){ // pbr

  Np=N(nanos);
  NMp=NM(nanos);
  Sp=S(nanos);
 
 for (int i=1;i<=nanos_proy;i++){

  Npplus=Np(nedades)*Sp(nedades);
  Np(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
  Np(nedades)+=Npplus;

  Np(1)= pR*mfexp(log_Rmed);//mean(Reclutas); // proyecta con R promedio
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

  if(mceval_phase())
    {
    ofstream out("joaq.doris",ios::out);
    //out << "kk" << 666 << "ll" << endl;
    out << CTP(1) << endl;
    out.close();
   }


REPORT_SECTION
  
  report << "Cruceros obs_pred" << endl;
  report << Bcrucero << endl;
  report << Bcrucero_pred << endl;

  report << "Capturas obs_pred" << endl;
  report << Desemb << endl;
  report << Desemb_pred << endl;

  report << "CPUE obs_pred" << endl;
  report << CPUE << endl;
  report << CPUE_pred << endl;

  report << "Btot " << endl;
  report << Btot << endl;

  report << "Bdesov " << endl;
  report << SD << endl;

  report << "Bvul " << endl;
  report << Bv << endl;
  
  report << "RPRs" << endl;
  report << RPR << endl;
  report << RPR2 << endl;
  report << "Reclutas" <<endl;
  report <<Reclutas<<endl;

  report << "N " << endl;
  report << N << endl;

  report << "F" << endl;
  report << F << endl;

  report << "Sflo" << endl;
  report << Sflo << endl;

  report << "p_obs " << endl;
  report << pobs << endl;
  report << "p_pred " << endl;
  report << ppred << endl;
  report << "pcruceros_obs " << endl;
  report << pobs_cru << endl;
  report << "pcruceros_pred " << endl;
  report << ppred_cru << endl;
  report << "Scruceros " << endl;
  report << Scru << endl;
  report << "Lf_obs_pred" << endl;
  report << Tallas*trans(pobs)<< endl;
  report << Tallas*trans(ppred)<< endl;
  report << "Lc_obs_pred" << endl;
  report << Tallas*trans(pobs_cru)<< endl;
  report << Tallas*trans(ppred_cru)<< endl;
  report << "_____________________________ " << endl;
  report << "log-like " << endl;
  report << likeval << endl;
  report << "_____________________________ " << endl;
  report << "Crecimiento " << endl;
  report <<"Loo    k    Lo "<<endl;
  report << Linf <<" "<< k <<"  "<< Lo <<endl;
  report << "Edades " << endl;
  report << edades << endl;
  report << "Talla_edad " << endl;
  report << mu_edad << endl;
  report << "$igma_edad " << endl;
  report << sigma_edad << endl;

  report << "_____________________________ " << endl;
  report << "Capturabilidad " << endl;
  report <<"Crucero"<<endl;
  report << qCru<<endl;
  report <<"BD virginal"<<endl;
  report << SDo2<<endl;
  report << "_____________________________ " << endl;
  report << "Matriz_p" << endl;
  report << Prob_talla << endl;
  report << "_____________________________ " << endl;
  report << "BD_proy"<<endl;
  report << SDp <<endl;      
  report << "BT_proy"<<endl;
  report << Bp << endl;
  report << "Y_proy"<<endl;
  report << YTP << endl;
  report << "Mortalidad natural" << endl;
  report << M <<endl;
  report << "Captura a la edad" << endl;
  report << Cedad <<endl;
  report << "RPR proyectado" << endl;
  report << RPRp << endl;
  report << SDp(nanos_proy)/SDo2 << endl;
  report << "FPBR" <<endl;
  report << Fpbr <<endl;

 
FUNCTION Eval_mcmc
  if(reporte_mcmc == 0)
  mcmc_report<<"CTP1 CTP2 CTP3 CTP4 CTP5 BDp1_fin BDp2_fin BDp3_fin BDp4_fin BDp5_fin"<<endl;
  mcmc_report<<YTP(1,1)<<" "<<YTP(1,2)<<" "<<YTP(1,3)<<" "<<YTP(1,4)<<" "<<YTP(1,5)<<
     " "<<SDp(nanos_proy,1)<<" "<<SDp(nanos_proy,2)<<" "<<SDp(nanos_proy,3)<<
     " "<<SDp(nanos_proy,4)<<" "<<SDp(nanos_proy,5)<<endl;

  reporte_mcmc++;


FINAL_SECTION

 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<"*********************************************"<<endl;

//GLOBALS_SECTION
//  #include  <admodel.h>
//  ofstream mcmc_report("mcmc.txt");
