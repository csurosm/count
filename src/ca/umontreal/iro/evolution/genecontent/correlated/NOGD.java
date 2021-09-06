
package ca.umontreal.iro.evolution.genecontent.correlated;


import ca.umontreal.iro.evolution.TreeNode;
import ca.umontreal.iro.evolution.genecontent.SimpleOccurrenceTable;

/**
 * Methods for computing non-orthologous gene displacement likelihoods
 * 
 * @author Louis Philippe B. Bouvrette and Mikl&oacute;s Cs&#369;r&ouml;s
 */
public class NOGD 
{
    
    /*aligne les indices c   */
    public static void mapping(TreeNode[] node){
        for (int i=0; i<node.length; i++){
             node[i].setId(i);
        }
    }
    
    /**
     * 
     * Carries out Felsenstein's peeling algorithm
     * 
     * @param leafProfile character value at leaves
     * @param node array of nodes in postorder traversal
     * @param subst substitution probabilities
     * @param root_distribution distribution of character state at root
     * @return likelihood of the profile
     */
    private static double getLikelihood(int[] leafProfile, TreeNode[] node, Exponentials[][] subst, double[] root_distribution)
    {
        double[][] likelihood =new double[node.length][8];
        double[][] P = new double[8][8];
        int leaf_idx=0;

        for (int node_idx =0; node_idx<node.length; node_idx++)
        {
            TreeNode N = node[node_idx];
            if (N.isLeaf())
                { // remplir likelihood[node_idx][] selon profile[leaf_idx]
                    for(int i=0; i<8; i++)
                    {
                        likelihood[node_idx][i] = 0;
                    }
                    likelihood[node_idx][leafProfile[leaf_idx]] = 1;
                    leaf_idx++;
                } else
                {// faire la rŽcurrence utilisant les likelihood[][] sur les enfants
                // L(u,a) = (Sum p(a->b) * L(v,b)) * (Sum p(a->b) * L(v',b)) * ...
                    int num_children = N.getNumChildren();
                    for(int a=0; a<8; a++){
                        likelihood[node_idx][a] = 1.0;
                        for (int i=0; i<num_children; i++){
                            TreeNode C = N.getChild(i);
                            double Lv = 0.0;
                            double t = C.getLength();
                            int cidx = C.getId();
                            for(int b=0; b<8; b++){
                               P[a][b] = subst[a][b].eval(t);
                               if (false && i==0){ // debug
                                System.out.println("#*NOGD.gCL "+node_idx+"/"+N.newickName()+"\ta "+a+"\tb "+b+"\tP "+P[a][b]+"\tt "+t);
                               }
                               Lv = Lv + P[a][b] * likelihood[cidx][b];
                            }
                            likelihood[node_idx][a] = likelihood[node_idx][a] * Lv;
                        }
                    }
                }
            if (false) { // debug
                System.out.print("#*NOGD.gCL "+node_idx+"/"+N.newickName());
                for (int a=0; a<8; a++)
                    System.out.print("\t"+Integer.toBinaryString(a+8).substring(1)+":"+likelihood[node_idx][a]);
                System.out.println();
            }
                        
        }

        int root_idx = node.length-1;
        double treecL = 0.0;
        for (int a=0; a<8; a++)
            treecL += likelihood[root_idx][a] * root_distribution[a];
    	return treecL;        
    }
    
    /**
     * Distribution of root state. X=1 with certainty, and Y-Y' are lost in tandem 
     * 
     * @param lambda gain rate
     * @param nu loss rate 
     * @param delta dispensability rate 
     * @return root state distribution
     */
    private static double[] getCorrelatedRootDistribution(double lambda, double nu, double delta)
    {
        double[] pi = new double[8];
        double omega = lambda*lambda+2.*lambda*delta+nu*delta;
        pi[0] = 0;
        pi[1] = 0;
        pi[2] = 0;
        pi[3] = 0;
        pi[4] = nu*delta/omega;
        pi[5] = delta*lambda/omega;
        pi[6] = delta*lambda/omega;
        pi[7] = lambda*lambda/omega;
        return pi;
    }
    
    /**
     * Distribution of root state. X=1 with certainty, and Y-Y' are lost independently
     * 
     * @param lambda gain rate
     * @param delta loss rate 
     * @return root state distribution
     */
    private static double[] getUncorrelatedRootDistribution(double lambda, double delta)
    {
        double[] pi = new double[8];
        double alpha = lambda + delta;
        double omega = alpha*alpha;
        pi[0] = 0;
        pi[1] = 0;
        pi[2] = 0;
        pi[3] = 0;
        pi[4] = delta*delta/omega;
        pi[5] = delta*lambda/omega;
        pi[6] = delta*lambda/omega;
        pi[7] = lambda*lambda/omega;
        return pi;
    }
    
    
    public static double getCorrelatedLikelihood(int[] leafProfile, TreeNode[] node, double mu, double lambda, double nu, double delta)
    {
        //System.out.println("#*NOGD.gCL mu "+mu+"\tlm "+lambda+"\tnu "+nu+"\tdl "+delta);
        Exponentials[][] S1 = getCorrelatedExponentials(lambda,nu,delta);
        Exponentials[][] S0 = getUncorrelatedExponentials(lambda,delta);
        Exponentials[][] Ct = getAllCorrelatedTransitions(S0, S1, mu); 
        double[] root_distribution = getCorrelatedRootDistribution(lambda, nu, delta);
        return getLikelihood(leafProfile, node, Ct, root_distribution);
    }
    
    public static double getUncorrelatedLikelihood(int[] leafProfile, TreeNode[] node, double mu, double lambda, double delta)
    {
        Exponentials[][] S0 = getUncorrelatedExponentials(lambda,delta);
        Exponentials[][] Ut = getAllUncorrelatedTransitions(S0, mu); 
        double[] root_distribution = getUncorrelatedRootDistribution(lambda, delta);
        return getLikelihood(leafProfile, node, Ut, root_distribution);
    }
    
    public static double getUnconditionalLikelihood(int[] leafProfile, TreeNode[] node, double mu, double lambda, double nu, double delta)
    {
        Exponentials[][] S1 = getCorrelatedExponentials(lambda,nu,delta);
        Exponentials[][] Ut = getAllUncorrelatedTransitions(S1, mu); 
        double[] root_distribution = getCorrelatedRootDistribution(lambda, nu, delta);
        return getLikelihood(leafProfile, node, Ut, root_distribution);
    }
    
    public static double getMixtureLikelihood( int[] leafProfile, TreeNode[] node, double mu, double lambda, double nu, double delta, double theta)
    {
        Exponentials[][] S1 = getCorrelatedExponentials(lambda,nu,delta);
        Exponentials[][] Ut = getAllUncorrelatedTransitions(S1, mu); 
        double[] root_distribution_unconditional = getCorrelatedRootDistribution(lambda, nu, delta);
        Exponentials[][] S0 = getUncorrelatedExponentials(lambda,delta);
        Exponentials[][] Ct = getAllCorrelatedTransitions(S0, S1, mu); 
        double[] root_distribution_conditional = getCorrelatedRootDistribution(lambda, nu, delta);
        
        Exponentials[][] Mt = new Exponentials[8][8];
        double[] root_distribution = new double[8];
        Exponentials Theta = new Exponentials(theta);
        Exponentials Theta_1 = new Exponentials(1.0-theta);
        for (int i=0; i<8; i++)
        {
            root_distribution[i] = theta*root_distribution_conditional[i]+(1.-theta)*root_distribution_unconditional[i];
            for (int j=0; j<8; j++)
            {
                Mt[i][j] = Ct[i][j].multiply(Theta).add(Ut[i][j].multiply(Theta_1));
            }
        }
        
        return getLikelihood(leafProfile, node, Mt, root_distribution);
    }
            
    
     
     
     /*Felsenstein peeling algorithm
      public static double getCorrelatedLikelihood(int[] leafProfile, TreeNode[] node, double mu, double lambda, double nu, double delta){
          //System.out.println("#*NOGD.gCL mu "+mu+"\tlm "+lambda+"\tnu "+nu+"\tdl "+delta);
        double[][] likelihood =new double[node.length][8];
        double[][] P = new double[8][8];
        int leaf_idx=0;
        Exponentials[][] S1 = getCorrelatedExponentials(lambda,nu,delta);
        Exponentials[][] S0 = getUncorrelatedExponentials(lambda,delta);
        Exponentials[][] Ct = getAllCorrelatedTransitions(S0, S1, mu); 
        for (int node_idx =0; node_idx<node.length; node_idx++)
        {
            TreeNode N = node[node_idx];
            if (N.isLeaf())
                { // remplir likelihood[node_idx][] selon profile[leaf_idx]
                    for(int i=0; i<8; i++){
                        likelihood[node_idx][i] = 0;
                    }
                    likelihood[node_idx][leafProfile[leaf_idx]] = 1;
                    leaf_idx++;
                }else
                {// faire la rŽcurrence utilisant les likelihood[][] sur les enfants
                // L(u,a) = (Sum p(a->b) * L(v,b)) * (Sum p(a->b) * L(v',b)) * ...
                    int num_children = N.getNumChildren();
                    for(int a=0; a<8; a++){
                        likelihood[node_idx][a] = 1.0;
                        for (int i=0; i<num_children; i++){
                            TreeNode C = N.getChild(i);
                            double Lv = 0.0;
                            double t = C.getLength();
                            int cidx = C.getId();
                            for(int b=0; b<8; b++){
                               P[a][b] = Ct[a][b].eval(t);
                               if (false && i==0){ // debug
                                System.out.println("#*NOGD.gCL "+node_idx+"/"+N.newickName()+"\ta "+a+"\tb "+b+"\tP "+P[a][b]+"\tt "+t);
                               }
                               Lv = Lv + P[a][b] * likelihood[cidx][b];
                            }
                            likelihood[node_idx][a] = likelihood[node_idx][a] * Lv;
                        }
                    }
                }
            if (false) { // debug
                System.out.print("#*NOGD.gCL "+node_idx+"/"+N.newickName());
                for (int a=0; a<8; a++)
                    System.out.print("\t"+Integer.toBinaryString(a+8).substring(1)+":"+likelihood[node_idx][a]);
                System.out.println();
            }
                        
        }
        double treecL = treeLikelihood(likelihood, node.length, lambda, nu, delta);
    	return treecL;
    }    

      public static double getUncorrelatedLikelihood(int[] leafProfile, TreeNode[] node, double mu, double lambda, double delta){
        double[][] likelihood =new double[node.length][8];
        double[][] P = new double[8][8];
        int leaf_idx=0;
        Exponentials[][] S0 = getUncorrelatedExponentials(lambda,delta);
        Exponentials[][] Ut = getAllUncorrelatedTransitions(S0, mu); 
        for (int node_idx =0; node_idx<node.length; node_idx++)
        {
            TreeNode N = node[node_idx];
            if (N.isLeaf())
                { // remplir likelihood[node_idx][] selon profile[leaf_idx]
                    for(int i=0; i<8; i++){
                        likelihood[node_idx][i] = 0;
                    }
                    likelihood[node_idx][leafProfile[leaf_idx]] = 1;
                    leaf_idx++;
                }else
                {// faire la rŽcurrence utilisant les likelihood[][] sur les enfants
                // L(u,a) = (Sum p(a->b) * L(v,b)) * (Sum p(a->b) * L(v',b)) * ...
                    int num_children = N.getNumChildren();
                    for(int a=0; a<8; a++){
                        likelihood[node_idx][a] = 1.0;
                        for (int i=0; i<num_children; i++){
                            TreeNode C = N.getChild(i);
                            double Lv = 0.0;
                            double t = C.getLength();
                            int cidx = C.getId();
                            for(int b=0; b<8; b++){
                               P[a][b] = Ut[a][b].eval(t);
                               Lv = Lv + P[a][b] * likelihood[cidx][b];
                            }
                            likelihood[node_idx][a] = likelihood[node_idx][a] * Lv;
                        }
                    }
                }
        }  
        double treeuL = treeLikelihood2(likelihood, node.length, lambda, delta);
    	return treeuL;
    }
      */

   /* calcule la vraisemblance de toutes les observations aux feuilles pour l'arbre correler
      public static double treeLikelihood(double[][] L, int root, double lambda, double nu, double delta){
          double treeL = 0.0;
          double[] pi = new double[8];
          double omega = lambda*lambda+2.*lambda*delta+nu*delta;
          pi[0] = 0;
          pi[1] = 0;
          pi[2] = 0;
          pi[3] = 0;
          pi[4] = nu*delta/omega;
          pi[5] = delta*lambda/omega;
          pi[6] = delta*lambda/omega;
          pi[7] = lambda*lambda/omega;
          for(int a = 0; a<8;a++){
          //    System.out.println(root-1);
            treeL = treeL + L[root-1][a] * pi[a];
           
          }
          return treeL; 
      }
    */
      /* calcule la vraisemblance de toutes les observations aux feuilles pour l'arbre non-correler
     public static double treeLikelihood2(double[][] L, int root, double lambda, double delta){
          double treeL = 0.0;
          double[] pi = new double[8];
          double alpha = lambda + delta;
          double omega = alpha*alpha;
          pi[0] = 0;
          pi[1] = 0;
          pi[2] = 0;
          pi[3] = 0;
          pi[4] = delta*delta/omega;
          pi[5] = delta*lambda/omega;
          pi[6] = delta*lambda/omega;
          pi[7] = lambda*lambda/omega;
          for(int a = 0; a<8;a++){
            treeL = treeL + L[root-1][a] * pi[a];
          }
          return treeL; 
      }*/

  public static Exponentials[][] getAllCorrelatedTransitions(Exponentials[][] S0, Exponentials[][] S1, double mu)
    {
        Exponentials[][] S = new Exponentials[8][8];
        // X transitions 0->0
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                S[i][j] = S0[i][j]; 

        // X transitions 0->1: constant 0
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                S[i][j+4] = new Exponentials(0.0);
            
        // X transitions 1->0
         for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
            {
                S[i+4][j] = new Exponentials(0.0);
                for (int k=0; k<4; k++)
                {
                   Exponentials Sk = S1[i][k].convolution(mu, S0[k][j]);
                   Exponentials Splus =  S[i+4][j].add(Sk);
                   //System.out.println("#*NOGD.gACT i "+i+"\tj "+j+"\tk "+k+"\tSk "+Sk);
                   S[i+4][j]=Splus;
                }
             }
        
        Exponentials loss_prob;
        {
            double[] cof = new double[1];
            cof[0] = 1.0;
            double[] exp = new double[1];
            exp[0] = mu;
            loss_prob = new Exponentials(cof,exp);
        }
        
       // X transitions 1->1
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
               S[i+4][j+4] = loss_prob.multiply(S1[i][j]);
                       
        return S;
    }
  
   public static Exponentials[][] getAllUncorrelatedTransitions(Exponentials[][] S0, double mu){
    Exponentials[][] S = new Exponentials[8][8];
               // X transitions 0->0
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                S[i][j] = S0[i][j]; 

        // X transitions 0->1: constant 0
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                S[i][j+4] = new Exponentials(0.0);
    
            Exponentials loss_prob;
        {
            double[] cof = new double[1];
            cof[0] = 1.0;
            double[] exp = new double[1];
            exp[0] = mu;
            loss_prob = new Exponentials(cof,exp);
        }
 
         Exponentials One = new Exponentials(1.0);
        // X transitions 1->0
         for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
            {
                for (int k=0; k<4; k++)
                    
                   S[i+4][j] = One.subtract(loss_prob).multiply(S0[i][j]);
             }
        
       
       // X transitions 1->1
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
               S[i+4][j+4] = loss_prob.multiply(S0[i][j]);
    return S;
   }
   
   public static Exponentials[][] getCorrelatedExponentials(double lambda, double nu, double delta)
    {
        double alpha = lambda+nu;
        double gamma = Math.sqrt(alpha*alpha+4.*(delta-nu)*(delta-lambda));
        double u = 3.*lambda+2.*delta+nu+gamma;
        double v = 3.*lambda+2.*delta+nu-gamma;
        double omega = lambda*lambda+2.*lambda*delta+nu*delta;
        double denom1 = (8.*omega*gamma);
	double denom2 = (2.*gamma*omega);
	double z = lambda / denom1;
        Exponentials pi00 = new Exponentials(nu*delta/omega);       //verified
        Exponentials pi01 = new Exponentials(delta*lambda/omega);   //verified 
        Exponentials pi10 = new Exponentials(delta*lambda/omega);   //verified
        Exponentials pi11 = new Exponentials(lambda*lambda/omega);  //verified
        Exponentials a,b,c,d,e,f,g;
	//a
        {
	    double[] cof = new double[2];
            cof[0] = (-3.*lambda+2.*delta-nu-gamma)*v*z;
            cof[1] = (3.*lambda-2.*delta+nu-gamma)*u*z;
            double[] exp = new double[2];
            exp[0] = u/2.0;
            exp[1] = v/2.0;
            a = new Exponentials(cof,exp);
        }
        //b
        {
            double y = lambda * lambda / denom2;
            double[] cof = new double[2];
            cof[0] = y*v;
            cof[1] = y*-u;
            double[] exp = new double[2];
            exp[0] = u/2.0;
            exp[1] = v/2.0;
            b = new Exponentials(cof,exp);
            

        }
	//c
        {
            double y = -(nu/ denom1);
            double[] cof = new double[2];
            cof[0] = (3.*lambda-2.*delta+nu+gamma)*v*y;
            cof[1] = (-3.*lambda+2.*delta-nu+gamma)*u*y;
            double[] exp = new double[2];
            exp[0] = u/2.0;
            exp[1] = v/2.0;
	    c = new Exponentials(cof,exp);
        }
	//d
        {
            double y = 1 / denom1;
	    double w = (delta - lambda) * (delta - lambda);
            double[] cof = new double[2];
            double gpa = ((gamma + alpha) / 2) * ((gamma + alpha) / 2);
            double gma = ((gamma - alpha) / 2) * ((gamma - alpha) / 2);
            cof[0] = -(w-gpa)*y*v;
            cof[1] = (w-gma)*y*u;
            double[] exp = new double[2];
            exp[0] = u/2.0;
            exp[1] = v/2.0;
                double[] fcof = new double[1];
                fcof[0] = -1.0/2.0;
                double[] fexp = new double[1];
                fexp[0] = alpha;
                Exponentials d_1 = new Exponentials(fcof,fexp);
	    Exponentials d_2 = new Exponentials(cof,exp);
            d = d_1.add(d_2);    
        }
	//e
	{
	    double[] cof = new double[2];
	    cof[0] = (lambda - 2.*delta - nu - gamma)*v*z;
            cof[1] = (-lambda + 2.*delta + nu - gamma)*u*z;
            double[] exp = new double[2];
            exp[0] = u/2.0;
            exp[1] = v/2.0;
            e = new Exponentials(cof,exp);
	}
	//f
       {
            double y = nu * delta / denom2;
            double[] cof = new double[2];
            cof[0] = y*v;
            cof[1] = y*-u;
            double[] exp = new double[2];
            exp[0] = u/2.0;
            exp[1] = v/2.0;
            f = new Exponentials(cof,exp);
        }
	//g
        {
	    double y = delta / denom1;
	    double[] cof = new double[2];
            cof[0] = (lambda - 2.*delta - nu - gamma)*v*y;
            cof[1] = (-lambda + 2.*delta + nu - gamma)*u*y;
            double[] exp = new double[2];
            exp[0] = u/2.0;
            exp[1] = v/2.0;
            g = new Exponentials(cof,exp);
        }
	//correlated transition matrix
        Exponentials[][] S1 = new Exponentials[4][4];
        S1[0][0] = pi00.subtract(a).subtract(a).subtract(b);//000->000
        S1[0][1] = pi01.add(a);//000->001
        S1[0][2] = pi10.add(a);//000->010
	S1[0][3] = pi11.add(b);//000->011
        S1[1][0] = pi00.add(c);//001->000
        S1[1][1] = pi01.subtract(c).subtract(d).subtract(e);//001->001
	S1[1][2] = pi10.add(d);//001->010
        S1[1][3] = pi11.add(e);//001->011
	S1[2][0] = pi00.add(c);//010->000
        S1[2][1] = pi01.add(d);//010->001
        S1[2][2] = pi10.subtract(c).subtract(d).subtract(e);//010->010
	S1[2][3] = pi11.add(e);//010->011
        S1[3][0] = pi00.add(f);//011->000
	S1[3][1] = pi01.add(g);//011->001
        S1[3][2] = pi10.add(g);//011->010
        S1[3][3] = pi11.subtract(f).subtract(g).subtract(g);//011->011
      /*
        double S00 = S1[0][0].eval(512);
        System.out.println(S00);
        double S01 = S1[0][1].eval(512);
        System.out.println(S01);
        double S02 = S1[0][2].eval(512);
        System.out.println(S02);
        double S03 = S1[0][3].eval(512);
        System.out.println(S03);
        double S10 = S1[1][0].eval(512);
        System.out.println(S10);
        double S11 = S1[1][1].eval(512);
        System.out.println(S11);
        double S12 = S1[1][2].eval(512);
        System.out.println(S12);
        double S13 = S1[1][3].eval(512);
        System.out.println(S13);
        double S20 = S1[2][0].eval(512);
        System.out.println(S20);
        double S21 = S1[2][1].eval(512);
        System.out.println(S21);
        double S22 = S1[2][2].eval(512);
        System.out.println(S22);
        double S23 = S1[2][3].eval(512);
        System.out.println(S23);
        double S30 = S1[3][0].eval(512);
        System.out.println(S30);
        double S31 = S1[3][1].eval(512);
        System.out.println(S31);
        double S32 = S1[3][2].eval(512);
        System.out.println(S32);
        double S33 = S1[3][3].eval(512);
        System.out.println(S33);
        System.out.println("-------------------"); 
       */
        return S1;
    }

       public static Exponentials[][] getUncorrelatedExponentials(double lambda, double delta)
    {
        double alpha = lambda+delta;
        double omega = alpha*alpha;
        Exponentials pi00 = new Exponentials(delta*delta/omega);
        Exponentials pi01 = new Exponentials(delta*lambda/omega);
        Exponentials pi10 = new Exponentials(delta*lambda/omega);
        Exponentials pi11 = new Exponentials(lambda*lambda/omega);
	Exponentials a,b,c,d,e,f,g;
	//a
	{
            double[] cof = new double[2];
            cof[0] = (lambda*lambda-lambda*delta)/omega;
            cof[1] = -lambda*lambda/omega;
            double[] exp = new double[2];
            exp[0] = alpha;
            exp[1] = 2.*alpha;
            a = new Exponentials(cof,exp);
        }
        //b
	{
            double[] cof = new double[2];
            cof[0] = (-2*(lambda*lambda))/omega;
            cof[1] = lambda*lambda/omega;
            double[] exp = new double[2];
            exp[0] = alpha;
            exp[1] = 2.*alpha;
            b = new Exponentials(cof,exp);
        }
	//c
	{
            double[] cof = new double[2];
            cof[0] = (lambda*delta-delta*delta)/omega;
            cof[1] = -lambda*delta/omega;
            double[] exp = new double[2];
            exp[0] = alpha;
            exp[1] = 2.*alpha;
            c = new Exponentials(cof,exp);
        }
	//d
	{
            double[] cof = new double[2];
            cof[0] = (-2*(lambda*delta))/omega;
            cof[1] = lambda*delta/omega;
            double[] exp = new double[2];
            exp[0] = alpha;
            exp[1] = 2.*alpha;
            d = new Exponentials(cof,exp);
        }
	//e
	{
            double[] cof = new double[2];
            cof[0] = (lambda*delta-lambda*lambda)/omega;
            cof[1] = -lambda*delta/omega;
            double[] exp = new double[2];
            exp[0] = alpha;
            exp[1] = 2.*alpha;
            e = new Exponentials(cof,exp);
        }
	//f
	{
            double[] cof = new double[2];
            cof[0] = (-2*(delta*delta))/omega;
            cof[1] = delta*delta/omega;
            double[] exp = new double[2];
            exp[0] = alpha;
            exp[1] = 2.*alpha;
            f = new Exponentials(cof,exp);
        }
	//g
	{
            double[] cof = new double[2];
            cof[0] = (delta*delta-delta*lambda)/omega;
            cof[1] = -delta*delta/omega;
            double[] exp = new double[2];
            exp[0] = alpha;
            exp[1] = 2.*alpha;
            g = new Exponentials(cof,exp);
        }
	//uncorrelated transition matrix
        Exponentials[][] S0 = new Exponentials[4][4];
        S0[0][0] = pi00.subtract(a).subtract(a).subtract(b);//000->000
	S0[0][1] = pi01.add(a);//000->001
        S0[0][2] = pi10.add(a);//000->010
	S0[0][3] = pi11.add(b);//000->011
        S0[1][0] = pi00.add(c);//001->000
        S0[1][1] = pi01.subtract(c).subtract(d).subtract(e);//001->001
	S0[1][2] = pi10.add(d);//001->010
        S0[1][3] = pi11.add(e);//001->011
	S0[2][0] = pi00.add(c);//010->000
        S0[2][1] = pi01.add(d);//010->001
        S0[2][2] = pi10.subtract(c).subtract(d).subtract(e);//010->010
	S0[2][3] = pi11.add(e);//010->011
        S0[3][0] = pi00.add(f);//011->000
	S0[3][1] = pi01.add(g);//011->001
        S0[3][2] = pi10.add(g);//011->010
        S0[3][3] = pi11.subtract(f).subtract(g).subtract(g);//011->011
  /*      
        double S00 = S0[0][0].eval(512);
        System.out.println(S00);
        double S01 = S0[0][1].eval(512);
        System.out.println(S01);
        double S02 = S0[0][2].eval(512);
        System.out.println(S02);
        double S03 = S0[0][3].eval(512);
        System.out.println(S03);
        double S10 = S0[1][0].eval(512);
        System.out.println(S10);
        double S11 = S0[1][1].eval(512);
        System.out.println(S11);
        double S12 = S0[1][2].eval(512);
        System.out.println(S12);
        double S13 = S0[1][3].eval(512);
        System.out.println(S13);
        double S20 = S0[2][0].eval(512);
        System.out.println(S20);
        double S21 = S0[2][1].eval(512);
        System.out.println(S21);
        double S22 = S0[2][2].eval(512);
        System.out.println(S22);
        double S23 = S0[2][3].eval(512);
        System.out.println(S23);
        double S30 = S0[3][0].eval(512);
        System.out.println(S30);
        double S31 = S0[3][1].eval(512);
        System.out.println(S31);
        double S32 = S0[3][2].eval(512);
        System.out.println(S32);
        double S33 = S0[3][3].eval(512);
        System.out.println(S33);
        */
        
        return S0;
    }      
      
    /**
     * Computes the transition matrix for the correlated evolution of three genes X, Y, Y'. The latter two
     * are redundant together when X is present. 
     * 
     * @author Mikl&oacute;s Cs&#369;r&ouml;s 
     * 
     * @param mu loss probability for gene X
     * @param lambda gain rate for genes Y,Y'
     * @param nu loss rate for X,Y,Y'
     * @param delta dispensability rate for Y,Y'
     * @return 8x8 matrix of transitions: row and column ordering by absence/presence pattern for XYY' 
     */
    public static Exponentials[][] getCorrelatedTransitions3(double mu, double lambda, double nu, double delta)
    {
        Exponentials[][] S1 = getCorrelatedTransitions2(lambda,nu,delta);
        Exponentials[][] S0 = getUncorrelatedTransitions2(lambda,delta);
        Exponentials[][] S = new Exponentials[8][8];
        // X transitions 0->0
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                S[i][j] = S0[i][j]; 
        
        // X transitions 0->1: constant 0
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                S[i][j+4] = new Exponentials(0.0); 
        
        // X transitions 1->0
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
            {
                S[i+4][j] = new Exponentials(0.0);
                for (int k=0; k<4; k++)
                    S[i+4][j] = S[i+4][j].add(S1[i][k].convolution(mu, S0[k][j]));
            }
        
        Exponentials loss_prob;
        {
            double[] cof = new double[1];
            cof[0] = 1.0;
            double[] exp = new double[1];
            exp[0] = mu;
            loss_prob = new Exponentials(cof,exp);
        }
        
        // X transitions 1->1
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                S[i+4][j+4] = loss_prob.multiply(S1[i][j]);
        
        
        return S;
    }

    /**
     * Computes the transition matrix for the correlated evolution of two genes YY'. 
     * 
     * @author Mikl&oacute;s Cs&#369;r&ouml;s 
     * 
     * @param lambda gain rate
     * @param nu loss rate 
     * @param delta dispensabilty rate
     * @return 4x4 matrix of transitions: row and column ordering by absence/presence pattern for YY'
     */
    public static Exponentials[][] getCorrelatedTransitions2(double lambda, double nu, double delta)
    {
        double alpha = lambda+nu;
        double gamma = Math.sqrt(alpha*alpha+4.*(delta-nu)*(delta-lambda));
        double u = 3.*lambda+2.*delta+nu+gamma;
        double v = 3.*lambda+2.*delta+nu-gamma;
        double omega = lambda*lambda+2.*lambda*delta+nu*delta;
        Exponentials pi00 = new Exponentials(nu*delta/omega);
        Exponentials pi01 = new Exponentials(delta*lambda/omega);
        Exponentials pi10 = new Exponentials(delta*lambda/omega);
        Exponentials pi11 = new Exponentials(lambda*lambda/omega);
        Exponentials a,b,c,d,e,f,g;
        {
            double z = lambda / (8.*omega*gamma);
            double[] cof = new double[2];
            cof[0] = (-3.*lambda+2.*delta-nu-gamma)*v*z;
            cof[1] = (3.*lambda-2.*delta+nu-gamma)*u*z;
            double[] exp = new double[2];
            exp[0] = u/2.0;
            exp[1] = v/2.0;
            a = new Exponentials(cof,exp);
        }
        // b,c,d,e,f,g
        Exponentials[][] S1 = new Exponentials[4][4];
        S1[0][1] = pi01.add(a);
        S1[0][2] = pi10.add(a);
        // ...
        return S1;
    }
    
    
    /**
     * Computes the transition matrix for the uncorrelated evolution of two genes YY'. 
     * 
     * @author Mikl&oacute;s Cs&#369;r&ouml;s 
     * 
     * @param lambda gain rate
     * @param delta loss rate
     * @return 4x4 matrix of transitions: row and column ordering by absence/presence pattern for YY'
     */
    public static Exponentials[][] getUncorrelatedTransitions2(double lambda, double delta)
    {
        double alpha = lambda+delta;
        double omega = alpha*alpha;
        Exponentials pi00 = new Exponentials(delta*delta/omega);
        Exponentials pi01 = new Exponentials(delta*lambda/omega);
        Exponentials pi10 = new Exponentials(delta*lambda/omega);
        Exponentials pi11 = new Exponentials(lambda*lambda/omega);
        
        Exponentials a,b,c,d,e,f,g;
        {
            double[] cof = new double[2];
            cof[0] = (lambda*lambda-lambda*delta)/omega;
            cof[1] = -lambda*lambda/omega;
            double[] exp = new double[2];
            exp[0] = alpha;
            exp[1] = 2.*alpha;
            a = new Exponentials(cof,exp);
        }
        // b,c,d,e,f,g
        Exponentials[][] S0 = new Exponentials[4][4];
        S0[0][1] = pi01.add(a);
        S0[0][2] = pi10.add(a);
        // ...
        
        return S0;
    }
    
    /**
     * Extracts the 0-7 profile for three families
     * 
     * @author Mikl&oacute;s Cs&#369;r&ouml;s 
     * 
     */
    public static int[] getTripleProfile(SimpleOccurrenceTable tbl, int x_idx, int y1_idx, int y2_idx)
    {
        int[] profile_x = tbl.getSizes(x_idx);
        int[] profile_y1 = tbl.getSizes(y1_idx);
        int[] profile_y2 = tbl.getSizes(y2_idx);
        
        int[] profile = new int[profile_x.length];
        for (int i=0; i<profile.length; i++)
        {
            profile[i]=0;
            if (profile_x[i]>0)
                profile[i]+=4;
            if (profile_y1[i]>0)
                profile[i]+=2;
            if (profile_y2[i]>0)
                profile[i]++;
        }
        return profile;
    }

}
