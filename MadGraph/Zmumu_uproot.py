import matplotlib.pyplot as plt
import os
import uproot
import awkward as ak

## path for ATLAS and CMS root files
folder = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/myExercises/Zmumu'
ATLAS_events = os.path.join(folder, 'Events/ATLAS/tag_1_delphes_events.root')
CMS_events = os.path.join(folder, 'Events/CMS/tag_1_delphes_events.root')
myATLAS_events = os.path.join(folder, 'Events/ATLAS/myATLAS.root')

def open_and_plot(file_root,label,output_dir):

    with uproot.open(f"{file_root}:Delphes") as events:
        
        # events.show() to check data type and names

        # choose the branches to save
        branches = events.arrays(["Muon.PT", "Muon.Eta", "Muon.Phi"], library='ak')

        muon_eta = ak.to_numpy(ak.flatten(branches['Muon.Eta']))
        muon_phi = ak.to_numpy(ak.flatten(branches['Muon.Phi']))


        fig, ax = plt.subplots()
        h = ax.hist2d(muon_eta, muon_phi, bins=[20,20])
        ax.set_title(rf'$pp \to Z \to \mu \mu$ {label}')  
        ax.set_xlabel(r'Muon $\eta$')
        ax.set_ylabel(r'Muon $\Phi$')
        cbar = fig.colorbar(h[3],ax=ax)
        cbar.set_label('Counts')

        plt.tight_layout()
        plot_path = os.path.join(output_dir, f"uproot_{label}.png")
        plt.savefig(plot_path, dpi=150)
        plt.close()

        return muon_eta, muon_phi         

###
if __name__ == '__main__':

    out_dir = '/home/mike/Physics/labParticles/MG5_aMC_v3_6_7/myScripts/plots'

    atlas_eta, atlas_phi = open_and_plot(ATLAS_events,'ATLAS',out_dir)
    cms_eta, cms_phi = open_and_plot(CMS_events,'CMS',out_dir)
    myatlas_eta, myatlas_phi = open_and_plot(myATLAS_events,'myATLAS',out_dir)

    plt.figure()
    plt.hist(atlas_eta,bins=30,label=r'ATLAS $\eta_{max}=2.5$',histtype='step')
    plt.hist(myatlas_eta,bins=30,label=r'myATLAS $\eta_{max}=2.0$',histtype='step')
    plt.title(r'$pp \to Z \to \mu \mu$ : ATLAS vs reduced track. eff.')
    plt.xlabel(r'Muon $\eta$')
    plt.ylabel('counts')
    plt.legend()
    #plt.show()
    plt.savefig(out_dir+'/uproot_myATLAS_diff.png',dpi=150)
    