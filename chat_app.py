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

# Initialize chat history if not exists
if "messages" not in st.session_state:
    st.session_state.messages = []

# Display chat history
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.write(message["content"])

# Chat input handler
prompt = st.chat_input("Enter SMILES or ask a question")

if prompt:
    # Add user message to chat
    st.chat_message("user").write(prompt)
    st.session_state.messages.append({"role": "user", "content": prompt})

    # Try to parse as SMILES first
    try:
        mol = Chem.MolFromSmiles(prompt)

        if mol is not None:
            # Its a valid SMILES string
            with st.chat_message("assistant"):
                st.write("Molecule parsed successfully!")

                # Display structure
                img = Draw.MolToImage(mol)
                st.image(img, width=400)

                # Calculate and display properties
                props = calculate_properties(mol)
                if props:
                    for name, value in props.items():
                        st.write(f"**{name}:** {value:.2f}")

                # Save assistant response
                st.session_state.messages.append({
                    "role": "assistant",
                    "content": "Molecular structure and properties displayed"
                })

        else:
            # Handle as a regular question
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
