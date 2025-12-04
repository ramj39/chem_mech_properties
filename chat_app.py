import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import py3Dmol

st.title("Molecular Property Calculator")

def calculate_properties(mol):
    """Calculate basic molecular properties"""
    try:
        return {
            "Molecular Weight": Descriptors.ExactMolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "H-Bond Donors": Descriptors.NumHDonors(mol),
            "H-Bond Acceptors": Descriptors.NumHAcceptors(mol),
            "Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
        }
    except Exception as e:
        st.error(f"Error calculating properties: {str(e)}")
        return None

def show_3d_structure(mol):
    """Render 3D structure using py3Dmol"""
    mb = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(mb, "mol")
    viewer.setStyle({"stick": {}})
    viewer.zoomTo()
    return viewer

# Initialize chat history
if "messages" not in st.session_state:
    st.session_state.messages = []

# Display chat history
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.write(message["content"])

# Chat input handler
prompt = st.chat_input("Enter SMILES or ask a question")

if prompt:
    st.chat_message("user").write(prompt)
    st.session_state.messages.append({"role": "user", "content": prompt})

    try:
        mol = Chem.MolFromSmiles(prompt)

        if mol is not None:
            with st.chat_message("assistant"):
                st.write("Molecule parsed successfully!")

                # 2D structure
                img = Draw.MolToImage(mol)
                st.image(img, width=400)

                # Properties
                props = calculate_properties(mol)
                if props:
                    for name, value in props.items():
                        st.write(f"**{name}:** {value:.2f}")

                # 3D structure
                viewer = show_3d_structure(mol)
                st.components.v1.html(viewer.render(), height=400)

                st.session_state.messages.append({
                    "role": "assistant",
                    "content": "Molecular structure (2D + 3D) and properties displayed"
                })

        else:
            with st.chat_message("assistant"):
                st.write("Invalid SMILES string. Please enter a valid SMILES or ask a chemistry question.")
                st.session_state.messages.append({
                    "role": "assistant",
                    "content": "Invalid SMILES string. Please enter a valid SMILES or ask a chemistry question."
                })

    except Exception as e:
        with st.chat_message("assistant"):
            st.error(f"Error processing input: {str(e)}")
            st.session_state.messages.append({
                "role": "assistant",
                "content": f"Error processing input: {str(e)}"
            })
