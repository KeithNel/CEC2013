/******************************************************************************
 * Version: 1.0
 * Last modified on: 30th July, 2014 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * c# port by Keith Nelson : keith_(AT)_cmind_(DOT)_org
 * ***************************************************************************/
namespace cec2013
{
    public interface IBoundsFn 
    {
	    int getDimension();
	    ClosedInterval.Double getBound( int dim );
    }
	
	public class ConstantBoundsFn : IBoundsFn
	{
		private int dim;
		private ClosedInterval.Double value;
		
		public ConstantBoundsFn( int dim, ClosedInterval.Double value ) {
			this.dim = dim;
			this.value = value;			
		}
				
		public int getDimension() {
			return dim;
		}

		public ClosedInterval.Double getBound(int dim) {
			return value;
		}
		
	}

	///////////////////////////////

	public class ExplicitBoundsFn : IBoundsFn
	{
		ClosedInterval.Double [] bounds;
		
		public ExplicitBoundsFn(params ClosedInterval.Double[] bounds ) {
			this.bounds = bounds;
		}
		
		public int getDimension() {
			return bounds.Length;
		}

		public ClosedInterval.Double getBound(int dim) {
			return bounds[ dim ];
		}
	}
}
