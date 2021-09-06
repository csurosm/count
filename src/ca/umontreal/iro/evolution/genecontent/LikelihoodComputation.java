
package ca.umontreal.iro.evolution.genecontent;

/**
 * A common interface used by likelihod computation algorithms
 * 
 * @author csuros
 */
import ca.umontreal.iro.matek.DiscreteDistribution;

public interface LikelihoodComputation 
{
    /**
     * Initializes the class with an occurrence table and a rate variation model.
     * 
     * @param table occurrence table
     * @param model rate variation model
     */
    public void init(OccurrenceTable table, RateVariation model);

    /**
     * The underlying occurrence table
     * 
     * @return the occurrence table
     */
    public OccurrenceTable getTable();
    
    /**
     * The underlying probablistic model
     * 
     * @return the rate variation model
     */
    public RateVariation getModel();
    

    /**
     * Recomputes the auxiliary values for a rate category
     * 
     * @param cidx category index in rate model
     */
    public void recomputeCategorySupport(int cidx);
    
    
    /**
     * Recomputes the auxiliary values for a rate category knowing that only one edge changed.
     * 
     * @param cidx category index in rate model
     * @param node_idx child end of the edge
     */
    public void recomputeCategorySupport(int cidx, int node_idx);
    

    /**
     * Recomputes the auxiliary values for profile-category pairs.
     * 
     * @param pidx profile index within the occurrence table
     * @param cidx category index in rate model
     */
    public void recomputeProfileSupport(int pidx, int cidx);
    
    /**
     * Recomputes the auxiliary values for profile-category pairs knowing that only one edge changed.
     * 
     * @param pidx profile index within the occurrence table
     * @param cidx category index in rate model
     * @param node_idx child end of the edge
     */
    public void recomputeProfileSupport(int pidx, int cidx, int node_idx);
    
    /**
     * Produces an array of conditional likelihoods at the root for a profile-rate category pair.
     * The returned array L[] is indexed by the number of surviving members at the root: L[n]
     * is the probablity of the population sizes at the leaves given that exactly n progenitors 
     * at the root have surviving descendants at at least one leaf. 
     * 
     * @param pidx profile index in occurrence table
     * @param cidx category index in rate variation model
     * @return array L[] of conditional likelihoods at the root
     */
    public double[] getRootLikelihoods(int pidx, int cidx);
    
    /**
     * Produces the prior probabilities for population at the root that can 
     * be multiplied with the root likelihoods to get the complete likelihood.
     * 
     * @param cidx category index in rate variation model
     * @return prior probabilities for root population sizes
     */
    public DiscreteDistribution getRootPrior(int cidx);
    
    /**
     * Produces the probability for phyletic profiles that are 
     * not included in the data set because they have too few members.
     * Data sets typically contain phyletic profiles only for those 
     * families that are represented in the lineages of interest. 
     * It may also be that only such profiles are included that 
     * are present in at least two terminal taxa: this is often the case 
     * in ortholog grouping.  
     * 
     * @param min_lineages number of terminal taxa at which a family must be present, may be 0, 1 or 2
     * @return probability of generating a profile with fewer lineages 
     */
    public double getAbsentProfileProbability(int min_lineages);


    public double[] getNodeLikelihoods(int node_idx, int pidx, int cidx);

    
}
