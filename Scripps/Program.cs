using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using pwiz.CLI.proteome;
using System.IO;
//using SourceList = System.Collections.Generic.Set<IDPicker.SourceInfo>;

namespace IDPicker
{
    // Summary:
    //     represents a peptide (sequence of amino acids)
    public class Peptide : IDisposable
    {
        public Peptide();
        public Peptide(string sequence);

        // Summary:
        //     returns the sequence of amino acids using standard single character symbols
        public string sequence { get; set; }
        public int uniqueCount { get; set; }
        public void Add(int uniqueCount, string peptideSeq)
        {
            this.uniqueCount = uniqueCount;
            this.sequence = peptideSeq;
        }

        
    }
    

    class Pep2Pro
    {
        

        /// <summary>
        /// find unique strings from a string list
        /// or simply, remove repeative strings from a string list
        /// </summary>
        /// <param name="inputList"></param>
        /// <returns></returns>
        public static List<string> removeDuplicate(List<string> inputList)
        {
            Dictionary<string, int> uniqueStore = new Dictionary<string, int>();
            List<string> finalList = new List<string>();

            foreach (var currValue in inputList)
            {
                if (!uniqueStore.ContainsKey(currValue))
                {
                    uniqueStore.Add(currValue, 0);
                    finalList.Add(currValue);
                }
            }
            return finalList;
        }

        

        /// <summary>
        /// This function removes protein clusters that do not have a 
        /// required number of additional unique peptides.
        /// </summary>
        /// <param name="minAdditionalPeptides"></param>
        public void filterByMinimumCoveringSet(int minAdditionalPeptides, ProteinList proteins)
        {

            if (minAdditionalPeptides == 0)
                return;

            if (statusOutput == null)
            {
                // Get the minimum protein set that can explain all the peptides
                List<ProteinInfo> proteinsToRemove = new List<ProteinInfo>();
                foreach (ProteinList.MapPair proItr in proteins)
                {
                    // Remove the protein cluster if it doesn't have the required number
                    // of unique peptides
                    if (proItr.Value.proteinGroup.uniquePeptideCount < minAdditionalPeptides)
                    {
                        proteinsToRemove.Add(proItr.Value);
                    }
                    else
                    {
                        proItr.Value.proteinGroup = null;
                        foreach (ResultInfo r in proItr.Value.results)
                            r.peptideGroup = null;
                    }
                }

                foreach (ProteinInfo pro in proteinsToRemove)
                    removeProtein(pro);
            }
            else
            {
                // Same as above with status reporting.
                List<ProteinInfo> proteinsToRemove = new List<ProteinInfo>();
                ProteinList.Enumerator itr = proteins.GetEnumerator(); itr.MoveNext();
                for (int i = 0; i < proteins.Count; ++i, itr.MoveNext())
                {
                    ProteinInfo pro = itr.Current.Value;
                    reportStatus("testing protein " + (i + 1) + " of " + proteins.Count, i + 1 == proteins.Count);

                    if (pro.proteinGroup.uniquePeptideCount < minAdditionalPeptides)
                        proteinsToRemove.Add(pro);
                    else
                    {
                        pro.proteinGroup = null;
                        foreach (ResultInfo r in pro.results)
                            r.peptideGroup = null;
                    }
                }

                for (int i = 0; i < proteinsToRemove.Count; ++i)
                {
                    reportStatus("removing protein " + (i + 1) + " of " + proteinsToRemove.Count, i + 1 == proteinsToRemove.Count);
                    removeProtein(proteinsToRemove[i]);
                }
            }

            assembleProteinGroups();
            assemblePeptideGroups();
            assembleClusters();

            // we lost the actual value by regrouping, but we know the minimum
            foreach (ProteinGroupInfo proGroup in proteinGroups)
                proGroup.uniquePeptideCount = minAdditionalPeptides;
        }

        /////////////////////////////////////////////////////////
        #region protein group, peptide group, and cluster assembly
        /// <summary>
        /// This functions adds proteins to groups. Each protein found
        /// across multiple samples are collapsed into groups based on
        /// whether they share same peptides or not. 
        /// </summary>
        public void assembleProteinGroups()
        {
            proteinGroups.Clear();
            if (statusOutput == null)
            {
                // Add each protein to the list of protein groups while 
                // remembering the protein group to which the protein belongs
                foreach (ProteinList.MapPair proItr in proteins)
                {
                    proItr.Value.proteinGroup = proteinGroups.addProtein(proItr.Value);
                }
            }
            else
            {
                // Do the same as above while reporting the status to the user.
                ProteinList.Enumerator itr = proteins.GetEnumerator(); itr.MoveNext();
                for (int i = 0; i < proteins.Count; ++i, itr.MoveNext())
                {
                    ProteinInfo pro = itr.Current.Value;
                    reportStatus("adding protein " + (i + 1) + " of " + proteins.Count, i + 1 == proteins.Count);
                    pro.proteinGroup = proteinGroups.addProtein(pro);
                }
            }
        }

        /// <summary>
        /// This function creates peptide groups. A peptide group contains
        /// peptide idenfications that matched to same proteins.
        /// </summary>
        public void assemblePeptideGroups()
        {
            peptideGroups.Clear();
            if (statusOutput == null)
            {
                // Take each result set and create peptide groups.
                // A peptide group contains peptide identifications
                // that matched to the same peptide across multiple
                // samples and multiple proteins.
                foreach (ResultInfo r in results)
                {
                    peptideGroups.addPeptide(r);
                }
            }
            else
            {
                // Do the same as above while reporting the status to the user
                ResultList.Enumerator itr = results.GetEnumerator(); itr.MoveNext();
                for (int i = 0; i < results.Count; ++i, itr.MoveNext())
                {
                    ResultInfo r = itr.Current;
                    reportStatus("adding result " + (i + 1) + " of " + results.Count, i + 1 == results.Count);
                    peptideGroups.addPeptide(r);
                }
            }
        }

        /// <summary>
        /// This function takes a protein group, assigns all the peptides and results to the supplied
        /// cluster. It then takes all the peptides and the results in the cluster, traces all the 
        /// proteins mapped to those peptides and the resutls and adds them to the same cluster.
        /// </summary>
        /// <param name="proGroup">A <see cref="IDPicker.ProteinGroupInfo"/> object containing
        /// protein identification results</param>
        /// <param name="c">A <see cref="IDPicker.ClusterInfo"/> object</param>
        /// 
        private void recursivelyAssignProteinGroupsToCluster(ProteinGroupInfo proGroup, ClusterInfo c)
        {
            reportStatus("assigning protein group " + proGroup.id + " to cluster " + c.id, true);

            if (proGroup.cluster > 0)
            {
                if (proGroup.cluster != c.id)
                    throw new InvalidDataException("protein groups that are connected are assigned to different clusters");

                return;
            }

            // Add the protein group to the cluster
            proGroup.cluster = c.id;
            c.proteinGroups.Add(proGroup);

            // For each protein in the group
            foreach (ProteinList.MapPair proItr in proGroup.proteins)
            {
                c.proteins[proItr.Value.locus] = proItr.Value;
                // Get the results and assign them to the cluster.
                // Also assign the corresponding peptides to the 
                // same cluster.
                foreach (ResultInfo r in proItr.Value.results)
                {
                    c.results.Add(r);
                    c.peptideGroups.Add(r.peptideGroup);
                    r.peptideGroup.cluster = c.id;
                }
            }

            // recursively add all "cousin" protein groups to the same cluster
            foreach (ResultInfo r in proGroup.results)
                foreach (ProteinGroupInfo cousinProGroup in r.peptideGroup.proteinGroups)
                {
                    if (!ReferenceEquals(cousinProGroup, proGroup) && cousinProGroup.cluster == 0)
                        recursivelyAssignProteinGroupsToCluster(cousinProGroup, c);
                    else if (cousinProGroup.cluster != c.id)
                        throw new InvalidDataException("protein groups that are connected are assigned to different clusters");
                }
        }

        /// <summary>
        /// This function takes peptide and protein groups and 
        /// assembles them together. A peptide group contains
        /// all peptide identifications that mapped to same 
        /// protein and a protein group contains all proteins that
        /// mapped to same peptides. The function generate clusters 
        /// and assign the proteins that are connected to the same 
        /// peptides to a single cluster. The proteins that map to
        /// a subset of the peptides in the generated cluster are 
        /// also assigned to the same cluster in a recursive fashion.
        /// </summary>
        public void assembleClusters()
        {
            clusters.Clear();

            if (statusOutput == null)
            {
                // Get each protein
                foreach (ProteinGroupInfo proGroup in proteinGroups)
                {
                    foreach (ResultInfo r in proGroup.results)
                    {
                        // Assign each peptide to corresponding group
                        proGroup.peptideGroups.Add(r.peptideGroup);
                        // Set the peptide to protein group mapping
                        r.peptideGroup.proteinGroups.Add(proGroup);
                    }
                }
            }
            else
            {
                // Do the same as above while reporting the status to the user
                ProteinGroupList.Enumerator itr = proteinGroups.GetEnumerator();
                for (int i = 0; i < proteinGroups.Count; ++i)
                {
                    itr.MoveNext();
                    ProteinGroupInfo proGroup = itr.Current;
                    IEnumerator<ResultInfo> itr2 = proGroup.results.GetEnumerator();
                    for (int j = 0; j < proGroup.results.Count; ++j)
                    {
                        itr2.MoveNext();
                        ResultInfo r = itr2.Current;
                        reportStatus("linking peptide group " + (j + 1) + " of " + proGroup.results.Count + " to protein group " + (i + 1) + " of " + proteinGroups.Count, j + 1 == proGroup.results.Count);
                        proGroup.peptideGroups.Add(r.peptideGroup);
                        r.peptideGroup.proteinGroups.Add(proGroup);
                    }
                }
            }

            // Set the protein and peptide group cluster identifications to 0.
            foreach (ProteinGroupInfo proGroup in proteinGroups)
            {
                proGroup.cluster = 0;
            }
            foreach (PeptideGroupInfo pepGroup in peptideGroups)
            {
                pepGroup.cluster = 0;
            }

            // For each protein group
            foreach (ProteinGroupInfo proGroup in proteinGroups)
            {
                if (proGroup.cluster == 0)
                {
                    // Generate a new cluster and assign an incremental ID.
                    ClusterInfo c = new ClusterInfo();
                    c.id = clusters.Count + 1;
                    // Add the protein to the cluster. Also add proteins that
                    // matched to the peptides of the given protein to the 
                    // cluster
                    recursivelyAssignProteinGroupsToCluster(proGroup, c);
                    clusters.Add(c);
                }
            }

            // Sort the clusters
            clusters.Sort(ClusterList.SortDescendingBySequencesThenSpectra);

            // Assign cluster idenfication numbers. These numbers signify
            // how many proteins in a cluster share peptides, and also 
            // how many protein are identified per peptide cluster.
            for (int i = 0; i < clusters.Count; ++i)
            {
                clusters[i].id = i + 1;
                foreach (ProteinGroupInfo proGroup in clusters[i].proteinGroups)
                    proGroup.cluster = i + 1;
                foreach (PeptideGroupInfo pepGroup in clusters[i].peptideGroups)
                    pepGroup.cluster = i + 1;
            }
        }
        public void assembleMinimumCoveringSet(ClusterInfo c)
        {
            if (c.proteinGroups.Count == 1) // degenerate case
            {
                foreach (ProteinGroupInfo proGroup in c.proteinGroups)
                    proGroup.uniquePeptideCount = int.MaxValue; // value is n/a
                return;
            }

            /*Set<ResultInfo> clusterResults = new Set<ResultInfo>( c.results );
            ProteinGroupList clusterGroups = new ProteinGroupList();
            foreach( ProteinGroupInfo proGroup in c.proteinGroups )
                clusterGroups.Add( proGroup );
            //Console.WriteLine(); 
            while( clusterResults.Count > 0 )
            {
                List<ProteinGroupInfo> minRemainingResults = new List<ProteinGroupInfo>();
                int minRemainingResultCount = clusterResults.Count;
                //int n = 0;
                //Console.WriteLine( "groups: " + clusterGroups.Count + "; results: " + clusterResults.Count );
                foreach( ProteinGroupInfo proGroup in clusterGroups )
                {
                    //Console.Write( n++ + " of " + clusterGroups.Count + "\r" );
                    int count = clusterResults.Count;
                    foreach( ResultInfo r in proGroup.results )
                        if( clusterResults.Contains( r ) )
                            --count;
                    if( count <= minRemainingResultCount )
                    {
                        if( count < minRemainingResultCount )
                            minRemainingResults.Clear();
                        minRemainingResults.Add( proGroup );
                    }
                }

                ProteinGroupInfo mostGreedyGroup = minRemainingResults[0];
                minRemainingResults.Clear();
                int oldCount = clusterResults.Count;
                clusterResults.Subtract( mostGreedyGroup.results );
                if( clusterResults.Count >= oldCount )
                {
                    Console.Error.WriteLine( "Something has gone terribly wrong!" );
                    System.Diagnostics.Process.GetCurrentProcess().Kill();
                }

                mostGreedyGroup.minSet = true;
                clusterGroups.Remove( mostGreedyGroup );
            }*/

            // Get the results in the cluster
            Set<ResultInfo> clusterResults = new Set<ResultInfo>(c.results);
            // Get the protein groups in the cluster
            ProteinGroupList clusterGroups = new ProteinGroupList();
            foreach (ProteinGroupInfo proGroup in c.proteinGroups)
                clusterGroups.Add(proGroup);
            //Console.WriteLine(); 
            // while there are results in the cluster
            while (clusterResults.Count > 0)
            {
                // Maps the number of results to a protein group
                Map<int, List<ProteinGroupInfo>> remainingResults = new Map<int, List<ProteinGroupInfo>>();
                //int n = 0;
                //Console.WriteLine( "groups: " + clusterGroups.Count + "; results: " + clusterResults.Count );
                // Iterate through protein groups
                foreach (ProteinGroupInfo proGroup in clusterGroups)
                {
                    //Console.Write( n++ + " of " + clusterGroups.Count + "\n" );
                    // Get the number of results in the cluster
                    int count = clusterResults.Count;
                    // Iterate over the cluster results and see how 
                    // many cluster group results can be explained
                    // by that protein group
                    foreach (ResultInfo r in proGroup.results)
                    {
                        if (clusterResults.Contains(r))
                            --count;
                    }
                    // Map the number of remaining results to that
                    // protein group
                    remainingResults[count].Add(proGroup);
                }

                // Take the first protein group that can explain the most results
                ProteinGroupInfo mostGreedyGroup = remainingResults.Values[0][0];
                // Subtract its results from the cluster results
                mostGreedyGroup.uniquePeptideCount = clusterResults.Count - remainingResults.Keys[0];
                clusterResults.Subtract(mostGreedyGroup.results);
                // Remove the most greedy group from the cluster groups
                clusterGroups.Remove(mostGreedyGroup);
            }
        }
        #endregion

        static void Main(string[] args)
        {
            Peptide peptideList = new Peptide();

            List<string> peptideSeqList = new List<string>();
            string file = "C:\\Users\\Dong\\Desktop\\test.pep.list";
            TextReader tr = new StreamReader(file);
            string line;
            while ((line = tr.ReadLine()) != null)
            {
                string peptide;
                if (line.Contains("."))
                {
                    string[] strs = line.Split('.');
                    peptide = strs[1];
                }
                else
                    peptide = line;
                peptideSeqList.Add(peptide);
            }

            //remove duplicates of the peptides
            //make it a unique list
            peptideSeqList = removeDuplicate(peptideSeqList);
            //generate the peptide list
            //by putting the peptide
            foreach (var pep in peptideSeqList)
            { 
                //something wrong here!!!
                peptideList.Add(1, "");
            }

            //Pep2Pro pp = new Pep2Pro();
            string database = "C:\\Users\\Dong\\Desktop\\yates.fasta";
            ProteomeDataFile foo = new ProteomeDataFile(database);
            pwiz.CLI.proteome.ProteinList proteinList = foo.proteinList;
            List<Protein> finalProteinList = new List<Protein>();
            
            
            for (int i = 0; i < proteinList.size(); i++)
            {
                Protein protein = proteinList.protein(i, true);
                if (!protein.id.Contains("Reverse"))
                {
                    string sequence = protein.sequence;
                    //Console.WriteLine(protein.index);
                    foreach (var pep in peptideSeqList)
                    {
                        if (sequence.Contains(pep) && !finalProteinList.Contains(protein))
                        {
                            finalProteinList.Add(protein);
                        }
                    }
                    
                }
            }

            //test
            //results: got a protein list which corresponding all peptides
            //the "finalProteinList" is ready for filters. 
            //Console.WriteLine(proteinList.size());
            //Console.WriteLine(finalProteinList.Count);


            //filter by minDistinctPeptides

            //filter by filterByDistinctPeptides

            //filterByMinimumCoveringSet

            //done it!
        }
    }
}
