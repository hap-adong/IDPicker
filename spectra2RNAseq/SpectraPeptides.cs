using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace Spectra2RNAseq
{
    public class SpectraPeptides
    {
        public int spectra = 0;
        public int peptides = 0;
        public List <string> spectralist = new List <string> ();
        public List <string> peptideslist = new List<string> ();
        public SpectraPeptides(string path)
        {
            string name = Path.GetFileNameWithoutExtension(path);
            IDPicker.Workspace workspace = new IDPicker.Workspace();
            Package.loadWorkspace(ref workspace, path);
            foreach (IDPicker.SourceGroupList.MapPair groupItr in workspace.groups)
                foreach (IDPicker.SourceInfo source in groupItr.Value.getSources(true))
                    foreach (IDPicker.SpectrumList.MapPair sItr in source.spectra)
                    {
                        IDPicker.ResultInstance ri = sItr.Value.results[1];
                        IDPicker.VariantInfo vi = ri.info.peptides.Min;
                        string rawPepSequence = vi.ToString();
                        string pepSequence = vi.peptide.sequence;
                        string index = sItr.Value.id.index.ToString();
                        string spectraID = name + "." + index;
                        peptideslist.Add(pepSequence);
                        spectralist.Add(spectraID);
                    }
            peptideslist = Package.removeDuplicate(peptideslist);
            spectralist = Package.removeDuplicate(spectralist);
            peptides = peptideslist.Count;
            spectra = spectralist.Count;
        }
    }
}
