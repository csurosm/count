package count.model;

import count.ds.ProfileTable;
import count.io.GeneralizedFileReader;
import count.io.NewickParser;
import count.io.RateVariationParser;
import count.io.TableParser;
import count.matek.Logarithms;

public class MixedRateLikelihood
{
	public MixedRateLikelihood(MixedRateModel mixed_model, ProfileTable table)
	{
		this.mixed_model = mixed_model;
		this.class_factories = new Likelihood[mixed_model.getNumClasses()];
		this.class_index = new int[class_factories.length];
		this.log_class_prob = new double[class_factories.length]; // precomputed
		this.table = table;
		this.num_classes = 0;
	}
	
	private final MixedRateModel mixed_model;
	private final ProfileTable table;
	
	private final Likelihood[] class_factories;
	private final int[] class_index;
	private final double[] log_class_prob;
	private int num_classes;
	
	private void computeParameters()
	{
		num_classes=0;
		for (int c = 0; c<class_factories.length; c++)
		{
			double pc = mixed_model.getClassProbability(c);
			if (pc != 0.0)
			{
				RateModel.GLD rates =mixed_model.getClassModel(c);
				Likelihood factory = new Likelihood(rates, table);
				class_factories[num_classes] = factory;
				class_index[num_classes] = c;
				log_class_prob[num_classes] = Math.log(pc);
				num_classes++;
			}
		}
	}
	
	public double getProfileLL(int family_idx)
	{
		if (num_classes==0) computeParameters();
		
		double[] terms = new double[num_classes];
		for (int i=0; i<num_classes; i++)
		{
			Likelihood factory = class_factories[i];
			Likelihood.Profile P = factory.getProfileLikelihood(family_idx);
			double pLL = P.getLogLikelihood();
			terms[i] = pLL+ log_class_prob[i]; 
			System.out.println("#*MRL.gPLL "+family_idx+"\ti "+i+"\tLL "+pLL+"\t"+Math.exp(pLL));
		}
		double LL = Logarithms.sum(terms, num_classes);
		return LL;
	}
	
	public double getEmptyLL()
	{
		if (num_classes==0) computeParameters();
		
		double[] terms = new double[num_classes];
		for (int i=0; i<num_classes; i++)
		{
			Likelihood factory = class_factories[i];
			terms[i] = factory.getEmptyLL() + log_class_prob[i];
		}
		double LL = Logarithms.sum(terms, num_classes);
		return LL;
	}
	
	public double getLL()
	{
		double LL = 0.0;
		for (int f=0; f<table.getFamilyCount(); f++)
		{
			LL += getProfileLL(f);
		}
		return LL;
	}
	
	public double getCorrectedLL()
	{
		double LL = getLL();
		int F = table.getFamilyCount();
		double L0 = getEmptyLL();
		// LL-F*log(1-exp(L0))
		double x = -Math.expm1(L0);
		LL -= F*Math.log(x); 
		return LL;
	}


	public static void main(String[] args) throws Exception
	{
        if (args.length!=3)
        {
            System.err.println("Call as $0 phylogeny table rates");
            System.exit(2008);
        }
        String tree_file = args[0];
        String table_file = args[1];
        String rates_file = args[2];
        count.ds.Phylogeny tree = NewickParser.readTree(
        		GeneralizedFileReader.guessReaderForInput(tree_file));
        count.ds.AnnotatedTable table = TableParser.readTable(tree.getLeafNames(),
        		GeneralizedFileReader.guessReaderForInput(table_file),true);
        RateVariation input_model = RateVariationParser.readRates(
        		GeneralizedFileReader.guessReaderForInput(rates_file)
        		, tree);
        
//        System.out.println(RateVariationParser.printRates(input_model));
        
        MixedRateLikelihood mixed_factory = new MixedRateLikelihood(input_model, table);
        
        double LL = mixed_factory.getLL();
        double L0 = mixed_factory.getEmptyLL();
        double p0 = Math.exp(L0);
        double corrL = LL-table.getFamilyCount()*Math.log1p(-p0);
        System.out.println("log-likelihood "+LL+"\tL0 "+L0+"\tp0 "+Math.exp(L0)+"\tcorr "+corrL);
		
	}


}
