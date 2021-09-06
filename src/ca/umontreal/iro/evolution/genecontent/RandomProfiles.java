package ca.umontreal.iro.evolution.genecontent;


import java.util.Random;
import ca.umontreal.iro.banality.BasicExecutable;

/**
 * Produces random profiles.
 *
 * @author csuros
 */
public class RandomProfiles extends BasicExecutable
{
    public RandomProfiles(RateVariation rates)
    {
        RND = new Random();
        setRates(rates);
    }

    private RateVariation rates;
    private Random RND;

    private double[] cumulative_class_prob;
    private int[] class_indexes;

    private void setRates(RateVariation rates)
    {
        this.rates = rates;
        if (rates !=null)
            init();
    }

    private void init()
    {
        int num_classes = rates.getNumClasses();
        int num_possible_classes = 0;
        for (int cidx=0; cidx<num_classes; cidx++)
            if (rates.getClassProbability(cidx)>0.0)
                num_possible_classes++;
        cumulative_class_prob = new double[num_possible_classes];
        class_indexes = new int[num_possible_classes];
        for (int cidx=0, j=0; cidx<num_classes;cidx++)
        {
            double p = rates.getClassProbability(cidx);
            if (p >0.0)
            {
                class_indexes[j] = cidx;
                cumulative_class_prob[j] = p;
                if (j>0)
                    cumulative_class_prob[j]+=cumulative_class_prob[j-1];
                j++;
            }
        }
        cumulative_class_prob[num_possible_classes-1] = 1.0; // ill be used as a sentinel, make sure there is no numerical error
//        for (int j=0; j<num_possible_classes; j++)
//        {
//            System.out.println("#*RP.init class "+j+" idx "+class_indexes[j]+" cum "+cumulative_class_prob[j]);
//        }
    }

    private int rndClassIndex()
    {
        double r = RND.nextDouble();
        for (int j=0; j<class_indexes.length; j++)
            if (r<=cumulative_class_prob[j])
                return class_indexes[j];
        // should never get here ...
        return -1;
    }

    private OccurrenceTable randomTable(int num_profiles)
    {
        TreeWithRates main_tree = rates.getMainTree();
        int num_classes = rates.getNumClasses();
        int num_nodes = main_tree.getNumNodes();
        double[] cumul_root_prior = new double[0];
        double[][][] cumul_duplication = new double[num_classes][][];
        double[][][] cumul_gain = new double[num_classes][][];
        int num_edges = main_tree.getNumEdges();
        for (int cidx=0; cidx<num_classes; cidx++)
            if (rates.getClassProbability(cidx)>0.0)
            {
                cumul_duplication[cidx]=new double[num_edges][0];
                cumul_gain[cidx] = new double[num_edges][0];
            }

        int[][] histories= new int[num_profiles][];
        int[][] profiles = new int[num_profiles][];
        int[] categories = new int[num_profiles];

        for (int pidx=0; pidx<num_profiles; )
        {
            int cidx = rndClassIndex();
            TreeWithRates class_tree = rates.getRateTree(cidx);
            int[] family_size = new int[num_nodes];
            boolean all_zero = true;
            for (int node_idx=num_nodes-1; node_idx>=0; node_idx--)
            {
                NodeWithRates N = class_tree.getNode(node_idx);
                if (N.isRoot()) // node_idx==num_nodes-1
                {
                    double r = RND.nextDouble();
                    int k = -1;
                    int i = 0;
                    do
                    {
                        for (; i<cumul_root_prior.length; i++)
                        {
//                            System.out.println("#*RP.rT "+pidx+"\troot\t-1\ti " +i+"\tr "+r+"\tcdf "+cumul_root_prior[i]);
                            if (r<=cumul_root_prior[i])
                            {
                                k = i;
                                break;
                            }
                        }
                        if (k==-1)
                        {
                            double[] cdf = rates.getRootPrior().getDistribution(2*cumul_root_prior.length+2);
                            for (int j=1; j<cdf.length; j++)
                                cdf[j] += cdf[j-1];
                            cumul_root_prior = cdf;
                        }
                    } while (k == -1);
//                    System.out.println("#*RP.rT "+pidx+"\troot " +k+"\tr "+r);
                    family_size[node_idx] = k;
                } else // non-root node
                {
                    int parent_idx = class_tree.getParentIndex(node_idx);
                    int parent_size = family_size[parent_idx];
                    int child_size = 0;
                    // duplications
                    for (int gidx=0; gidx<parent_size; gidx++)
                    {
                        double r = RND.nextDouble();
                        int k = -1;
                        int i=0;
                        double[] cdf = cumul_duplication[cidx][node_idx];
                        do
                        {
                            int clen = cdf.length;
                            for (; i<clen; i++)
                                if (r<=cdf[i])
                                {
                                    k = i;
                                    break;
                                }
                            if (k==-1)
                            {
                                cdf = N.getDuplicationDistribution().getDistribution(2*clen+2);
                                for (int j=1; j<cdf.length; j++)
                                    cdf[j] += cdf[j-1];
                                cumul_duplication[cidx][node_idx] = cdf;
                            }
                        } while (k == -1);
                        child_size += k;
                    }
                    // gains
                    {
                        double r = RND.nextDouble();
                        int k = -1;
                        int i=0;
                        double[] cdf = cumul_gain[cidx][node_idx];
                        do
                        {
                            int clen = cdf.length;
                            for (; i<clen; i++)
                                if (r<=cdf[i])
                                {
                                    k = i;
                                    break;
                                }
                            if (k==-1)
                            {
                                cdf = N.getTransferDistribution().getDistribution(2*clen+2);
                                for (int j=1; j<cdf.length; j++)
                                    cdf[j] += cdf[j-1];
                                cumul_gain[cidx][node_idx] = cdf;
                            }
                        } while (k == -1);
                        child_size += k;
                    }
                    family_size[node_idx] = child_size;
                    if (child_size > 0 && N.isLeaf())
                        all_zero = false;
                } // if non-root node
            } // for all nodes
            if (!all_zero)
            {
                int num_leaves = main_tree.getNumLeaves();

                int[] pattern = new int[num_leaves];
                int leaf_idx = 0;
                for (int node_idx=0; node_idx<num_nodes; node_idx++)
                {
                    NodeWithRates N = main_tree.getNode(node_idx);
                    if (N.isLeaf())
                    {
                        pattern[leaf_idx] = family_size[node_idx];
                        leaf_idx++;
                    }
                }
                profiles[pidx] = pattern;
                histories[pidx] = family_size;
                categories[pidx] = cidx;
                pidx++;
            } // not all-0
        } // for all profiles

        OccurrenceTable tbl = new OccurrenceTable(main_tree.getLeaves());
        tbl.setTable(profiles);
        int[] prop_idx = new int[num_nodes];
        int cat_prop_idx = tbl.registerProperty("*Rate category");
        for (int node_idx=0; node_idx<num_nodes; node_idx++)
        {
            NodeWithRates N = main_tree.getNode(node_idx);
            if (!N.isLeaf())
                prop_idx[node_idx] = tbl.registerProperty("*"+N.getTaxonName());
        }
        String[] category_names = new String[num_classes];
        for (int cidx=0; cidx<num_classes; cidx++)
        {
            StringBuffer sb = new StringBuffer();
            sb.append(Integer.toString(cidx));
            sb.append(':');
            sb.append("E");
            sb.append(Integer.toString(rates.getEdgeLengthCategory(cidx)));
            sb.append(",G");
            sb.append(Integer.toString(rates.getTransferRateCategory(cidx)));
            sb.append(",L");
            sb.append(Integer.toString(rates.getLossRateCategory(cidx)));
            sb.append(",D");
            sb.append(Integer.toString(rates.getDuplicationRateCategory(cidx)));
            category_names[cidx] = sb.toString();
        }

        for (int pidx=0; pidx<num_profiles; pidx++)
        {
            int cidx = categories[pidx];
            tbl.setFamilyProperty(pidx, cat_prop_idx, category_names[cidx]);
            for (int node_idx=0; node_idx<num_nodes; node_idx++)
            {
                NodeWithRates N = main_tree.getNode(node_idx);
                if (!N.isLeaf())
                    tbl.setFamilyProperty(pidx, prop_idx[node_idx], Integer.toString(histories[pidx][node_idx]));

            }
        }

        return tbl;
    }

    private void go(String[] args) throws Exception
    {
        if (args.length != 3)
        {
            die(new IllegalArgumentException("Call as "+getClass().getCanonicalName()+" tree rates n"));
        }
        reportLaunch(args);
        String tree_file = args[0];
        String rates_file = args[1];
        int n = Integer.parseInt(args[2]);

        TreeWithRates input_tree = new TreeWithRates(ca.umontreal.iro.evolution.Parser.readNewick(new java.io.FileReader(tree_file)));
        RateVariation input_model = RateVariation.read(new java.io.FileReader(rates_file),input_tree);
        setRates(input_model);
        OccurrenceTable tbl = randomTable(n);
        System.out.println(tbl.getFormattedTable(true));


    }

    /**
     * Call as $0 tree rates num_profiles
     * @param args
     */
    public static void main(String[] args) throws Exception
    {
        RandomProfiles RP = new RandomProfiles(null);
        RP.go(args);
    }

}
