/******************************************************************************
 * Version: 1.0
 * Last modified on: 30th July, 2014 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * c# port by Keith Nelson : keith_(AT)_cmind_(DOT)_org
 * ***************************************************************************/
using System;

namespace cec2013
{
    public abstract class ClosedInterval : IComparable<ClosedInterval>, ICloneable
    {
        public ClosedInterval() {}
    
        public static bool isClosedInterval( double lower, double upper )
        {
	        return lower <= upper;
        }

        ///////////////////////////////
        
        public bool isEmpty()
        {
            return getLength() == 0;
        }

        ///////////////////////////////
    
        public abstract double getLower();
        public abstract double getUpper();

        ///////////////////////////////
            
        public bool contains( double x )
        {
	        return getLower() <= x && x <= getUpper();
        }

        public bool overlaps( ClosedInterval other )
        {
	        return contains( other.getLower() ) || contains( other.getUpper() );
        }

        public bool contains( ClosedInterval other )
        {
	        return contains( other.getLower() ) && contains( other.getUpper() );
        }

        public double getLength() 
        {
            double upper = getUpper();
            double lower = getLower();        
            return upper > lower ? upper - lower : 0;
        }

        public override int GetHashCode()
        {
            return (getLower().GetHashCode() >> 13) ^ getUpper().GetHashCode();
        }
	
	    public bool equals( Object o )
        {
            return o is ClosedInterval
                    && getLower() == ((ClosedInterval)o).getLower()
                        && getUpper() == ((ClosedInterval)o).getUpper();
        }
    
        public override string ToString()
        {
            return "[ " + getLower() + ", " + getUpper() + " ]";
        }
            
        public int CompareTo( ClosedInterval rhs )
        {
            if( getLower() < rhs.getLower() )
                return -1;
            else if( getLower() == rhs.getLower() )
            {
                if( getUpper() < rhs.getUpper() )
                    return -1;
                else if( getUpper() == rhs.getUpper() )
                    return 0;
                else
            	    return 1;
            }
            else
        	    return 1;        
        }    
    
        public abstract Object Clone();
    
        ///////////////////////////////
    
        public Double 
        intersection( ClosedInterval rhs )
        {
            Double result = new Double();

            if( overlaps( rhs ) )
            {
                double lower = getLower() > rhs.getLower() ? getLower() : rhs.getLower();
                double upper = getUpper() < rhs.getUpper() ? getUpper() : rhs.getUpper();
                result = new Double( lower, upper );
            }    

            return result;
        }
    
        ///////////////////////////////    

        public class Double : ClosedInterval
        {
            public double lower = 0;
            public double upper = -1;

            ///////////////////////////////

            public Double() { }

            public Double(double fl, double cl)
            {
                if (!isClosedInterval(fl, cl))
                    throw new ArgumentException("closed interval expected, found [" + fl + ',' + cl + ']');

                lower = fl;
                upper = cl;
            }

            public Double(Double rhs)
            {
                lower = rhs.lower;
                upper = rhs.upper;
            }

            public override double getLower() { return lower; }
            public override double getUpper() { return upper; }

            public new bool equals(Object o)
            {
                if (!(o is Double))
                    return false;

                Double rhs = (Double)o;
                return lower == rhs.lower && upper == rhs.upper;
            }

            public override object Clone()
            {
                return new Double(lower, upper);
            }
        }

        ///////////////////////////////    

        /*public class Int : ClosedInterval
        {
            public int lower = 0;
            public int  upper = -1;
        
            ///////////////////////////////
            
            public Int() {}

            public Int( int fl, int cl )
            {
                if (!isClosedInterval(fl, cl))
                    throw new ArgumentException("closed interval expected, found [" + fl + ',' + cl + ']');

                lower = fl;
                upper = cl;
            }

            public Int( Int rhs )
            {
                lower = rhs.lower;
                upper = rhs.upper;
            }

            public override double getLower() { return lower; }
            public override double getUpper() { return upper; }

            public new bool equals( Object o )
            {
                if( !( o is Int ) ) //Long ) )
                    return false;
        	
                //ClosedInterval.Long rhs = (ClosedInterval.Long)o;
                Int rhs = (Int)o;
                return lower == rhs.lower && upper == rhs.upper;
            }

            public override object Clone()
            {
                return new Int(lower, upper);
            }
        }
    
        ///////////////////////////////    

        public class Long : ClosedInterval
        {
            public long lower = 0;
            public long upper = -1;
        
            ///////////////////////////////
            
            public Long() {}

            public Long( long fl, long cl )
            {
                if (!isClosedInterval(fl, cl))
                    throw new ArgumentException("closed interval expected, found [" + fl + ',' + cl + ']');

                lower = fl;
                upper = cl;
            }

            public Long( Long rhs )
            {
                lower = rhs.lower;
                upper = rhs.upper;
            }

            public override double getLower() { return lower; }
            public override double getUpper() { return upper; }

            public new bool equals( Object o )
            {
                if( !( o is Long ) )
                    return false;

                Long rhs = (Long)o;
                return lower == rhs.lower && upper == rhs.upper;
            }

            public override object Clone()
            {
                return new Long(lower, upper);
            }
        }
        */
        ///////////////////////////////
    }
}
