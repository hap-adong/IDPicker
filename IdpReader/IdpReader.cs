using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Runtime.InteropServices;
using System.Xml;
using System.Text.RegularExpressions;
using System.Net;
using System.Reflection;
using System.ComponentModel;
using pwiz.CLI.analysis;
using pwiz.CLI.data;
using pwiz.CLI.msdata;
using pwiz.CLI.proteome;
//using SourceList = System.Collections.Generic.Set<IDPicker.SourceInfo>;



namespace IDPReader
{
    class IdpReader
    {

        public static void Main(string[] args)
        {
            string idpxml = "K:\\project vandy\\dataset\\ordinal regression 12-10-10\\human_DLD1_LTQ_nm.idpXML";
            string mzml = "K:\\project vandy\\dataset\\ordinal regression 12-10-10\\human_DLD1_LTQ.mzML";

            List<string> list = Package.idpReader(idpxml, mzml, 0.98, 2);

            TextWriter file = new StreamWriter(output);
            foreach (string ss in list)
            {
                file.WriteLine(ss);
            }



            /////////////////////////////////////////////////////
            //original methodology for idpReader
            /////////////////////////////////////////////////////
            //List<string> pepList = new List<string>();
            //string input = "Z:\\home\\dwang\\fragmentation\\RNA-Seq\\RKO\\UniversalPep.csv";
            //TextReader reader = new StreamReader(input);
            //string line;
            //while ((line = reader.ReadLine()) != null)
            //{
            //    pepList.Add(line);
            //}

            //string path = "Z:\\home\\dwang\\fragmentation\\RNA-Seq\\RKO\\myrimatch\\mam_121007n_RKO_200ug_IEF_";
            //string output = "Z:\\home\\dwang\\fragmentation\\RNA-Seq\\RKO\\myrimatch\\output.csv";
            //List<string> list = new List<string>();
            //for (int i = 1; i <= 3; i++)
            //{
            //    for (int j = 1; j <= 10; j++)
            //    {
            //        string let = "";
            //        string num = "";
            //        switch (i)
            //        {
            //            case 1:
            //                let = "A";
            //                break;
            //            case 2:
            //                let = "B";
            //                break;
            //            case 3:
            //                let = "C";
            //                break;
            //            default:
            //                Console.WriteLine("error");
            //                break;
            //        }
            //        switch (j)
            //        {
            //            case 1:
            //                num = "01";
            //                break;
            //            case 2:
            //                num = "02";
            //                break;
            //            case 3:
            //                num = "03";
            //                break;
            //            case 4:
            //                num = "04";
            //                break;
            //            case 5:
            //                num = "05";
            //                break;
            //            case 6:
            //                num = "06";
            //                break;
            //            case 7:
            //                num = "07";
            //                break;
            //            case 8:
            //                num = "08";
            //                break;
            //            case 9:
            //                num = "09";
            //                break;
            //            case 10:
            //                num = "10";
            //                break;
            //            default:
            //                Console.WriteLine("error");
            //                break;
            //        }
            //        string ext = let + num;
            //        string xml = path + ext + ".idpXML";
            //        string mzML = path + ext + ".mzML";
            //        list = Package.idpReader(xml, mzML, 0.98, 3, pepList, list);
            //        Console.WriteLine("one file done: " + (i * 10 + j));
            //    }
            //}//end for

            //TextWriter file = new StreamWriter(output);
            //file.WriteLine("nativeID,pepSequence,basePeak,TIC,b1,b2,b3,y1,y2,y3,len,bond,NBasic,CBasic,NR,NK,NH,NL,CR,CK,CH,CL,R,K,H,L,ambiguityLabel");
            //foreach (string ss in list)
            //{
            //    file.WriteLine(ss);
            //}
        }//end of main

    }//end IdpReader
}
