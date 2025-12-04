import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import py3Dmol

def calculate_properties(mol):
    """Calculate basic molecular properties"""
    return {
        "Molecular Weight": Descriptors.ExactMolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "H-Bond Donors": Descriptors.NumHDonors(mol),
        "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
        "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        #"heavy_atom_count":Descriptors.heavy_atom_count(mol),
    }

# Chat input handler
if prompt := st.chat_input("Enter SMILES or ask a question"):
    # Try to parse as SMILES first
    mol = Chem.MolFromSmiles(prompt)

    if mol is not None:
        # It's a valid SMILES string
        with st.chat_message("assistant"):
            st.write("Molecule parsed successfully!")

            # Display structure
            st.image(Draw.MolToImage(mol))

            # Calculate and display properties
            props = calculate_properties(mol)
            for name, value in props.items():
                st.write(f"**{name}:** {value:.2f}")

    else:
        # Handle as a regular question
        with st.chat_message("assistant"):
            st.write("Please ask your chemistry question...")
