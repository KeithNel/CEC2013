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
    public abstract class Func {

	    private readonly List<ClosedInterval.Double> bounds;
	
	    public Func( IBoundsFn boundsFn ) {
		    ClosedInterval.Double[] l = new ClosedInterval.Double[boundsFn.getDimension()]; 
		    for( int i=0; i<boundsFn.getDimension(); ++i )
			    l[i] = boundsFn.getBound( i );
		
		    bounds = new List<ClosedInterval.Double>( l );
	    }
	
	    public double evaluate( double [] x ) {
		    return evaluateCommon( x, true );
	    }
	
	    public double evaluate( List<double> x ) {
		    double [] xx = new double[x.Count];
            x.CopyTo(xx);
		    return evaluateCommon( xx, true );
	    }

	    private double evaluateCommon( double [] x, bool checkBounds ) {
		    if( checkBounds && !isInBounds( x ) )
			    throw new ArgumentException();

		    return doEvaluate( x );
	    }
	
	    public double cfuncEvaluate( double [] x ) {

		    return evaluateCommon( x, false );
	    }

        public virtual int getDimension() {
		    return bounds.Count;
	    }

	    public bool isInBounds( double [] x ) {
		    if( x.Length != getDimension() )
			    return false;
		
		    for( int i=0; i<x.Length; ++i )
			    if( !bounds[i].contains( x[ i ] ) )
				    return false;
		    return true;
	    }

        public List<ClosedInterval.Double> getBounds() {
            return bounds;
        }

        public abstract double doEvaluate(double[] x);
    }
}
