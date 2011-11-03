using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using pwiz.CLI.analysis;
using pwiz.CLI.data;
using pwiz.CLI.msdata;
using pwiz.CLI.proteome;

namespace IonRetrieval
{
    class NativeIDComparison
    {
        public static List<string> findCommon(List<string> A, List<string> B)
        {
            List<string> commonList = new List<string>();
            foreach (string a in A)
            {
                foreach (string b in B)
                {
                    if (a == b)
                    {
                        commonList.Add(a);
                    }
                }
            }

            return commonList;
        }

        public static List<string> IDList(string idpXML, string mzML)
        {
            List<string> list = new List<string>();
            MSDataFile foo = new MSDataFile(mzML);
            SpectrumList sl = foo.run.spectrumList;

            IDPicker.Workspace workspace = new IDPicker.Workspace();
            Package.loadWorkspace(ref workspace, idpXML);

            foreach (IDPicker.SourceGroupList.MapPair groupItr in workspace.groups)
                foreach (IDPicker.SourceInfo source in groupItr.Value.getSources(true))
                    foreach (IDPicker.SpectrumList.MapPair sItr in source.spectra)
                    {
                        IDPicker.ResultInstance ri = sItr.Value.results[1];
                        IDPicker.VariantInfo vi = ri.info.peptides.Min;
                        bool boolCharge = sItr.Value.id.charge.Equals(3);
                        if (boolCharge)
                        {
                            string rawPepSequence = vi.ToString();
                            string interpretation = vi.ToSimpleString();
                            // Look up the index with nativeID
                            object idOrIndex = null;
                            if (sItr.Value.nativeID != null && sItr.Value.nativeID.Length > 0)
                            {
                                idOrIndex = sItr.Value.nativeID;
                                list.Add(idOrIndex.ToString());
                            }
                        }//end if (boolcharge)
                    }//end foreach
            return list;
        }

        public static void nativeIDComp(string idpXML_nm, string idpXML_pm, string mzML)
        {
            List<string> nmNativeIDList = new List<string>();
            List<string> pmNativeIDList = new List<string>();
            List<string> nmExclusiveIDList = new List<string>();
            List<string> pmExclusiveIDList = new List<string>();

            nmNativeIDList = IDList(idpXML_nm, mzML);
            pmNativeIDList = IDList(idpXML_pm, mzML);
           
            List<string> commonList = findCommon(nmNativeIDList, pmNativeIDList);
            Console.WriteLine("common: " + commonList.Count);
            Console.WriteLine("nm: " + nmNativeIDList.Count);
            Console.WriteLine("pm: " + pmNativeIDList.Count);
            
            foreach (string nativeID in nmNativeIDList)
            {
                if (!commonList.Contains(nativeID)) nmExclusiveIDList.Add(nativeID);
            }

           
            foreach (string nativeID in pmNativeIDList)
            {
                if (!commonList.Contains(nativeID)) pmExclusiveIDList.Add(nativeID);
            }

            Console.WriteLine("MM model");
            foreach (string nativeID in nmExclusiveIDList) Console.WriteLine(nativeID);

            Console.WriteLine("Basophile model");
            foreach (string nativeID in pmExclusiveIDList) Console.WriteLine(nativeID);
        }
    }
}
