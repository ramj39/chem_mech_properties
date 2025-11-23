import streamlit as st
from mendeleev import element
st.markdown(
    """
    <style>
    body, .stApp {
        background: linear-gradient(45deg, #ff9a9e 0%, #fad0c4 99%,#fad0c4 100%);
        min-height: 100vh;
        background-attachment: fixed;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.title("🔍 Element QueryBot")

# User input
query = st.text_input("Enter element name or symbol (e.g., Fe or Iron):")

if query:
    try:
        el = element(query.strip().capitalize())

        st.subheader(f"Properties of {el.name} ({el.symbol})")
        st.markdown(f"- **Atomic Number**: {el.atomic_number}")
        st.markdown(f"- **Atomic Weight**: {el.atomic_weight}")
        st.markdown(f"- **Electron Configuration**: `{el.econf}`")
        st.markdown(f"- **Oxidation States**: {el.oxidation_states}")
        st.markdown(f"- **Electronegativity (Pauling)**: {el.en_pauling}")
        st.markdown(f"- **Ionization Energy**: {el.ionenergies.get(1)} eV")
        st.markdown(f"- **Melting Point**: {el.melting_point} K")
        st.markdown(f"- **Boiling Point**: {el.boiling_point} K")
        st.markdown(f"- **Density**: {el.density} g/cm³")
        st.markdown(f"- **Van der Waals Radius**: {el.vdw_radius} pm")
        st.markdown(f"- **Covalent Radius**: {el.covalent_radius} pm")
        st.markdown(f"- **Block**: {el.block}")
        st.markdown(f"- **Group**: {el.group_id}")
        st.markdown(f"- **Period**: {el.period}")
        st.markdown(f"- **Discovery Year**: {el.discovery_year}")
        st.markdown(f"- **Discovered By**: {el.discoverers}")

    except Exception as e:
        st.error(f"Could not find data for '{query}'. Try a valid element name or symbol.")
