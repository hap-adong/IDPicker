using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Data;
using System.IO;
namespace UPS
{
    class Program
    {
        static void Main(string[] args)
        {
            //===========================================================
            string folder = "Z:\\home\\dwang\\fragmentation\\UPS";
            string file_naive = folder + "\\mvh-naive\\klc_031308p_cptac_study6_6_QC1.idpXML";
            string file_baso = "Z:\\home\\dwang\\fragmentation\\UPS\\mvh-baso\\klc_031308p_cptac_study6_6_QC1.idpXML";
            List<string> peplist_naive = Package.PepSecurity(file_naive, 3);
            List<string> peplist_baso = Package.PepSecurity(file_baso, 3);
            peplist_naive = Package.removeDuplicate(peplist_naive);
            peplist_baso = Package.removeDuplicate(peplist_baso);

            int sup_naive = 0;
            int sup_baso = 0;
            string fasta = "Z:\\home\\dwang\\fragmentation\\UPS\\Sigma49-reverse.fasta";
            List<string> proteinList = Package.readDatabase(fasta);
            foreach (var peptide in peplist_naive)
            {
                foreach (var protein in proteinList)
                { 
                    if (protein.Contains(peptide))
                    {
                        sup_naive++;
                        break;
                    }
                }
            }

            foreach (var peptide in peplist_baso)
            {
                foreach (var protein in proteinList)
                {
                    if (protein.Contains(peptide))
                    {
                        sup_baso++;
                        break;
                    }
                }
            }

            Console.WriteLine(sup_naive);
            Console.WriteLine(sup_baso);

            //how is the overlap
            int overlap = 0;
            foreach (var pep_naive in peplist_naive)
            {
                foreach (var pep_baso in peplist_baso)
                {
                    if (pep_naive == pep_baso)
                    {
                        overlap++;
                        break;
                    }
                }
            }

            Console.WriteLine("overlap: " + overlap);
        }
    }
}
