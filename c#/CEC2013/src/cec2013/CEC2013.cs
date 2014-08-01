/******************************************************************************
 * Version: 1.0
 * Last modified on: 30th July, 2014 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * c# port by Keith Nelson : keith_(AT)_cmind_(DOT)_org
 * ***************************************************************************/
using System;
using System.Collections.Generic;

namespace cec2013
{
    /**
     * This class provides all necessary information and methods 
     * for the function set in the CEC2013 Competition on 
     * Niching methods for multimodal Optimisation 
     * 
     * @author Michael G. Epitropakis
     * c# port by Keith Nelson
     */
    public class CEC2013 {

	    /*Array of the final functions for the competition */
	    private static readonly List< Func > funcs_ = new List< Func >() {
			    new Funcs.FiveUnevenPeakTrap( new ConstantBoundsFn( 1, new ClosedInterval.Double(0.0, 30.0) ) ),
			    new Funcs.EqualMaxima( new ConstantBoundsFn( 1, new ClosedInterval.Double(0.0, 1.0) ) ),
			    new Funcs.UnevenDecreasingMaxima( new ConstantBoundsFn( 1, new ClosedInterval.Double(0.0, 1.0) ) ),
			    new Funcs.Himmelblau( new ConstantBoundsFn( 2, new ClosedInterval.Double( -6.0, 6.0 ) ) ),
			    new Funcs.SixHumpCamelBack( new ExplicitBoundsFn( 
					    new ClosedInterval.Double( -1.9, 1.9 ), 
					    new ClosedInterval.Double( -1.1, 1.1 ) ) ), 
				new Funcs.Shubert( new ConstantBoundsFn( 2, new ClosedInterval.Double( -10.0, 10.0 ) ) ),
				new Funcs.Vincent( new ConstantBoundsFn( 2, new ClosedInterval.Double( 0.25, 10.0 ) ) ),
				new Funcs.Shubert( new ConstantBoundsFn( 3, new ClosedInterval.Double( -10.0, 10.0 ) ) ),
				new Funcs.Vincent( new ConstantBoundsFn( 3, new ClosedInterval.Double( 0.25, 10.0 ) ) ),
				new Funcs.ModifiedRastriginAll( new ConstantBoundsFn( 2, new ClosedInterval.Double( 0.0, 1.0 ) ) ), 
				new Funcs.CF1( new ConstantBoundsFn( 2, new ClosedInterval.Double( -5, 5 ) ) ),
				new Funcs.CF2( new ConstantBoundsFn( 2, new ClosedInterval.Double( -5, 5 ) ) ),
				new Funcs.CF3( new ConstantBoundsFn( 2, new ClosedInterval.Double( -5, 5 ) ) ),
				new Funcs.CF3( new ConstantBoundsFn( 3, new ClosedInterval.Double( -5, 5 ) ) ),
				new Funcs.CF4( new ConstantBoundsFn( 3, new ClosedInterval.Double( -5, 5 ) ) ),
				new Funcs.CF3( new ConstantBoundsFn( 5, new ClosedInterval.Double( -5, 5 ) ) ),
				new Funcs.CF4( new ConstantBoundsFn( 5, new ClosedInterval.Double( -5, 5 ) ) ),
				new Funcs.CF3( new ConstantBoundsFn( 10, new ClosedInterval.Double( -5, 5 ) ) ),
				new Funcs.CF4( new ConstantBoundsFn( 10, new ClosedInterval.Double( -5, 5 ) ) ),
				new Funcs.CF4( new ConstantBoundsFn( 20, new ClosedInterval.Double( -5, 5 ) ) )
        };

	    /**
	     * Store the fitness of the global optima per function.
	     */
	    private static double [] fopt_ = new double[] { 
		    200.0, 1.0, 1.0, 200.0, 1.03163, 186.731, 1.0, 2709.0935, 1.0, -2.0, 
		    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	    /**
	     * Store the $\rho$ values per function.
	     */
	    private static double [] rho_ = new double[] { 
		    0.01, 0.01, 0.01, 0.01, 0.5, 0.5, 0.2, 0.5, 0.2, 0.01, 
		    0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };

	    /**
	     * Store the number of global optima per function. 
	     */
	    private static int [] nopt_ = new int[] { 
		    2, 5, 1, 4, 2, 18, 36, 81, 216, 12, 
		    6, 8, 6, 6, 8, 6, 8, 6, 8, 8};

	    /**
	     * Store the maximum number of function evaluation per function. 
	     */
	    private static long [] maxfes_ = new long [] { 
		    50000, 50000, 50000, 50000, 50000, 200000, 200000, 400000, 400000, 200000, 
		    200000, 200000, 200000, 400000, 400000, 400000, 400000, 400000, 400000, 400000 };

	    /**
	     * The default constructor for this class.
	     */
	    public CEC2013() { 	}

	    /**
	     * An easy way to get all the functions for the competition.
	     * @return a List of Func objects which corresponds to the 
	     * 			function set of the competition
	     */
	    public List< Func > getFunctions() { return funcs_; }

	    /**
	     * Get the fitness value of the global optima per function
	     * @param nfunc The number of function to be optimised (1..20) 
	     * @return the fitness value of the global optima of the nfunc function (1..20)
	     */
	    public static double getFitnessGoptima(int nfunc) { return fopt_[nfunc-1]; }

	    /**
	     * Get the dimension of the nfunc function (1..20)
	     * @param nfunc The number of function to be optimised (1..20)
	     * @return the dimension of the nfunc function (1..20)
	     */
	    public static int getDimension(int nfunc) { return funcs_[nfunc-1].getDimension(); }

	    /**
	     * Get the number of the global optima for the nfunc-th function
	     * @param nfunc The number of function to be optimised (1..20)
	     * @return the number of the global optima for the nfunc-th function
	     */
	    public static int getNoGoptima(int nfunc)  { return nopt_[nfunc-1]; }

	    /**
	     * Get the $\rho$ value for the nfunc-th function
	     * @param nfunc The number of function to be optimised (1..20)
	     * @return the $\rho$ value for the nfunc-th function
	     */
	    public static double getRho(int nfunc) { return rho_[nfunc-1]; }

	    /**
	     * Get the maximum number of function evaluations per function to be optimised
	     * @param nfunc The number of function to be optimised (1..20)
	     * @return the maximum number of function evaluations for the nfunc-th function
	     */
	    public static long getMaxfes(int nfunc) { return maxfes_[nfunc-1]; }

	    /**
	     * Calculate the seeds of a population (based on Algorithm 1 of the technical report)
	     * @param pop the population of the potential solutions an List of Lists
	     * 			e.g. population with size NPxDim: An List with size NP that contains 
	     * 				 Arraylists of "Dim" Doubles 
	     * @param radius the $\rho$ value
	     * @return the seeds of a population
	     */
	    private static List<List<double>> getSeeds(List<List<double>> pop, double radius)
	    {
		    List<List<double>> seeds = new List<List<double>>();
		    /* Determine the species seeds: 
		     * iterate through sorted population */
		    foreach(List<double> indi in pop) {
			    bool found = false;
			    /* Iterate seeds */
			    foreach (List<double> sindi in seeds) {
				    /* Calculate distance from seeds */
				    double sum = 0;
				    for (int j=0; j < sindi.Count; j++)
					    sum += ( indi[j]-sindi[j] ) * ( indi[j]-sindi[j] );

				    double dist = Math.Sqrt(sum);

				    /* If the Euclidean distance is less than the radius */
				    if (dist <= radius) {
					    found = true;
					    break;
				    }
			    } //seeds

			    /* If it is not similar to any other seed, then it is a new seed */
			    if (!found) { seeds.Add(indi); }
		    } //pop

		    return seeds;
	    }

	    public class CompareByFitness : IComparer<int> {
		    private List<double> fitness;

		    public CompareByFitness(List<double> fitness) {
			    this.fitness = fitness;
		    }

		    public int Compare(int o1, int o2) {
			    double f1 = fitness[o1];
			    double f2 = fitness[o2];

			    if (f1 < f2)
				    return 1;
			    else if (f1 > f2)
				    return -1;
			    else 
				    return 0;
		    }
	    }


	    /**
	     * This function calculates the number of global optima that exists inside a population  
	     * @param pop the population of the potential solutions
	     * @param seeds the calculated seeds of the population
	     * @param findex the index of the function that the user optimises 
	     * 			BE CAREFUL: the index should be in the range of (1...20)  
	     * @param accuracy the accuracy level for the experiment
	     * @param radius the radius $\rho$
	     * @return the number of global optima found within the pre-specified accuracy 
	     */
	    public static int howManyGlobalOptimaInPopulation(List<List<double>> pop,
			    List<List<double>> seeds, int findex,
			    double accuracy, double radius)
	    {
		    Func pFunc = CEC2013.funcs_[findex-1];
		    List<List<double>> sorted_pop = new List<List<double>>();

		    /* Evaluate pop */
		    List<double> fits = new List<double>();
		    foreach (List<double> indi in pop) 
			    fits.Add(pFunc.evaluate(indi));

		    /* sorted_pop: sort population by fitness */
		    CompareByFitness compare = new CompareByFitness(fits);
		    List<int> idx = new List<int>();

		    for (int i=0; i < pop.Count; i++) {
			    idx.Add(i);
		    }

            idx.Sort(compare);

		    foreach (int i in idx) {
			    sorted_pop.Add(pop[i]);
		    }

		    /* find seeds in the temp population */
		    List<List<double>> seedstmp = new List<List<double>>();
		    seedstmp.AddRange(getSeeds(sorted_pop, radius));

		    /* Based on the accuracy: check which seeds are global optimizers */
		    int count = 0;
		    foreach (List<double> indi in seedstmp) {
			    /* evaluate seed */
			    double seedFitness = pFunc.evaluate(indi);

			    /* |F_seed - F_goptimum| <= accuracy */
			    if ( Math.Abs(seedFitness - CEC2013.getFitnessGoptima(findex)) <= accuracy) {
				    seeds.Add(indi);
				    ++count;
			    }
			    /* save time */
			    if (count == CEC2013.getNoGoptima(findex)) break;
		    }

		    return count;
	    }
    }

// End ///////////////////////////////////////////////////////////////
}
