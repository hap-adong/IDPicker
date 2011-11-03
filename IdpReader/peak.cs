using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace IDPReader
{
    public class Peak : IComparable<Peak>
    {

        public double mz;
        public double rankOrIntensity;
<<<<<<< HEAD
        //added for the purpose of fragment ions for orbiorbi
        public int fragmentCharge;
=======
>>>>>>> 80fc9e47b4dafd28bcc96f478ae26543036b9ff6

        public Peak(double mass, double rnkOrIntens)
        {
            mz = mass;
            rankOrIntensity = rnkOrIntens;
        }

<<<<<<< HEAD
        public Peak(double mass, double rnkOrIntens, int z)
        {
            mz = mass;
            rankOrIntensity = rnkOrIntens;
            fragmentCharge = z;
        }

=======
>>>>>>> 80fc9e47b4dafd28bcc96f478ae26543036b9ff6
        public Peak(Peak pk)
        {
            mz = pk.mz;
            rankOrIntensity = pk.rankOrIntensity;
        }

        public Peak(double mass)
        {
            mz = mass;
            rankOrIntensity = -1;
        }

        public int CompareTo(Peak rhs)
        {
            if (mz < rhs.mz)
                return -1;
            else if (mz > rhs.mz)
                return 1;
            return 0;
        }
    }
}
