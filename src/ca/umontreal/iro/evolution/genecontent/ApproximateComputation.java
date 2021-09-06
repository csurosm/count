package ca.umontreal.iro.evolution.genecontent;

import ca.umontreal.iro.banality.Verbose;

/**
 *
 * @author csuros
 */
public class ApproximateComputation extends StableComputation
{
    public ApproximateComputation(OccurrenceTable table, RateVariation model)
    {
        super(table, model);
    }

    public ApproximateComputation(){ super(); }

    public static final double DEFAULT_CUTOFF_FACTOR = 3.0;
    public static final int DEFAULT_CUTOFF_THRESHOLD = 6;

    private double cutoff_factor = DEFAULT_CUTOFF_FACTOR; // this many times max. occurrence within subtree
    private int cutoff_threshold = DEFAULT_CUTOFF_THRESHOLD; // but at least this many

    public void setCutoffFactor(double factor)
    {
        this.cutoff_factor = factor;
    }

    public void setCutoffThreshold(int threshold)
    {
        this.cutoff_threshold = threshold;
    }

    @Override
    public void init(OccurrenceTable table, RateVariation model)
    {
        super.init(table, model);
    }

    @Override
    protected LowerConditionals newLowerConditionals(PhyleticProfile profile, Probabilities probs)
    {
        return new ApproximateLowerConditionals(profile, probs);
    }

    private class ApproximateLowerConditionals extends LowerConditionals
    {
        private ApproximateLowerConditionals(PhyleticProfile profile, Probabilities survival)
        {
            super(profile, survival);
        }

        @Override
        protected void initDataStructures()
        {
            // super.initDataStructures();
            NodeWithRates[] nodes = rate_tree.getDFT();
            edge_likelihood = new double[nodes.length-1][];
            node_likelihood = new double[nodes.length][];

            int max_occ = 0;
            for (int leaf_idx=0; leaf_idx<rate_tree.getNumLeaves(); leaf_idx++)
            {
                int occ = profile.get(leaf_idx);
                max_occ = (occ>max_occ?occ:max_occ);
            }
            int cutoff = Math.max(cutoff_threshold, (int)(cutoff_factor*max_occ));

            for (int node_idx=0; node_idx<nodes.length; node_idx++)
            {
                NodeWithRates N = rate_tree.getNode(node_idx);
                if (N.isLeaf())
                {
                    int n = profile.get(node_idx);
                    node_likelihood[node_idx] = new double[n+1];
                } else
                {
                    int num_children = N.getNumChildren();
                    int sum_m = 0;
                    for (int ci=0; ci<num_children; ci++)
                    {
                        int child_position = getChildForStep(node_idx,ci);
                        int child_idx = rate_tree.getChildIndex(node_idx, child_position);
                        int m = node_likelihood[child_idx].length-1;

                        sum_m += m;
                        edge_likelihood[child_idx] = new double[1+Math.min(sum_m,cutoff)];
                    }
                    int last_child_idx = rate_tree.getChildIndex(node_idx, getChildForStep(node_idx,num_children-1));
                    node_likelihood[node_idx] = edge_likelihood[last_child_idx];
                }
            }
            //if (Verbose.isVerbose())
            //{
            //    int[] s = profile.computeSubtreeSizes(rate_tree);
            //    for (int node_idx=0; node_idx<nodes.length; node_idx++)
            //    {
            //        Verbose.message("SC.LC.iDS node "+node_idx+"/"+rate_tree.getNode(node_idx).newickName()
            //                +"\tsize "+s[node_idx]+"\tnl[] "+node_likelihood[node_idx].length);
            //        for (int ci=0; ci<rate_tree.getNode(node_idx).getNumChildren(); ci++)
            //        {
            //            int child_idx = rate_tree.getChildIndex(node_idx, getChildForStep(node_idx,ci));
            //            Verbose.message("SC.LC.iDS edge "+ci+"/"+rate_tree.getNode(child_idx).newickName()
            //                    +"\tsize "+s[child_idx]+"\tel[] "+edge_likelihood[child_idx].length);
            //        }
            //    }

            //}

        }

    }

    private void go(String[] args) throws Exception
    {
        reportLaunch(args);
        reportOtherArguments("Maximum number of paralogs: "+ML.MAX_PARALOGS);
        reportOtherArguments("Minimum number of families: "+ML.MIN_PRESENT_LINEAGES);
        reportOtherArguments("Cutoff: "+this.cutoff_threshold+"\t| "+this.cutoff_factor+" * max occurrence");
        ML.testComputation(this, args);
    }

    public static void main(String[] args) throws Exception
    {
        Verbose.setVerbose(false);

        ApproximateComputation O = new ApproximateComputation();

        int num_switches = 0;
        try {
            while (args.length>2*num_switches && args[2*num_switches].startsWith("-"))
            {
                String arg_switch = args[2*num_switches].substring(1);
                if (arg_switch.equals("h"))
                    O.go(new String[0]); // will throw an Exception
                if (args.length==2*num_switches+1)
                    throw new IllegalArgumentException("Missing argument for switch "+args[2*num_switches]);
                String arg_value = args[2*num_switches+1];
                if (arg_switch.equals("v"))
                {
                    Verbose.setVerbose(arg_value.equals("true"));
                } else if (arg_switch.equals("max_paralogs"))
                {
                    ML.MAX_PARALOGS = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("min_lineages"))
                {
                    ML.MIN_PRESENT_LINEAGES = Integer.parseInt(arg_value);
                } else if (arg_switch.equals("factor"))
                {
                    O.cutoff_factor = Double.parseDouble(arg_value);
                } else if (arg_switch.equals("cutoff"))
                {
                    O.cutoff_threshold = Integer.parseInt(arg_value);
                }
                else
                    throw new IllegalArgumentException("Switch not recognized: '"+args[2*num_switches]+"'");

                num_switches++;
            }

            String[] rest=new String[args.length-2*num_switches];
            for (int j=0; j<rest.length; j++)
                rest[j]=args[2*num_switches+j];
            long T0 = System.currentTimeMillis();
            O.go(rest);
            long delta = System.currentTimeMillis()-T0;
            System.out.println("| Time: "+delta+" ms");
        } catch (Exception E)
        {
            die(E);
        }
    }

}
