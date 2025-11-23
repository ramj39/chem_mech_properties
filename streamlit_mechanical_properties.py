import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from math import pi, sqrt
st.markdown(
    """
    <style>
    body, .stApp {
        background: linear-gradient(135deg, #b3ffab 0%, #12fff7 100%);
        min-height: 100vh;
        background-attachment: fixed;
    }
    </style>
    """,
    unsafe_allow_html=True
)

# Set page configuration
st.set_page_config(
    page_title="Metallurgical Properties Calculator",
    page_icon="⚙",
    layout="wide"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .property-card {
        background-color: #f0f2f6;
        padding: 20px;
        border-radius: 10px;
        margin: 10px 0;
        border-left: 4px solid #1f77b4;
    }
    .result-box {
        background-color: #e8f4fd;
        padding: 15px;
        border-radius: 8px;
        margin: 10px 0;
        border: 1px solid #1f77b4;
    }
</style>
""", unsafe_allow_html=True)

# Material database with comprehensive properties
MATERIAL_DATABASE = {
    "AISI 1020 (Low Carbon Steel)": {
        "type": "Carbon Steel",
        "yield_strength": 350,
        "tensile_strength": 420,
        "elongation": 25,
        "hardness_hv": 130,
        "youngs_modulus": 200,
        "density": 7.87,
        "thermal_conductivity": 51.9,
        "specific_heat": 486,
        "melting_point": 1520,
        "coefficient_thermal_expansion": 11.7,
        "electrical_resistivity": 15.9,
        "poissons_ratio": 0.29,
        "carbon_content": 0.20
    },
    "AISI 1045 (Medium Carbon Steel)": {
        "type": "Carbon Steel",
        "yield_strength": 530,
        "tensile_strength": 625,
        "elongation": 12,
        "hardness_hv": 200,
        "youngs_modulus": 200,
        "density": 7.85,
        "thermal_conductivity": 50.2,
        "specific_heat": 486,
        "melting_point": 1515,
        "coefficient_thermal_expansion": 11.3,
        "electrical_resistivity": 17.3,
        "poissons_ratio": 0.29,
        "carbon_content": 0.45
    },
    "AISI 4140 (Alloy Steel)": {
        "type": "Alloy Steel",
        "yield_strength": 655,
        "tensile_strength": 1020,
        "elongation": 17,
        "hardness_hv": 300,
        "youngs_modulus": 205,
        "density": 7.85,
        "thermal_conductivity": 42.7,
        "specific_heat": 477,
        "melting_point": 1420,
        "coefficient_thermal_expansion": 12.3,
        "electrical_resistivity": 22.4,
        "poissons_ratio": 0.29,
        "carbon_content": 0.40
    },
    "AISI 304 (Stainless Steel)": {
        "type": "Stainless Steel",
        "yield_strength": 215,
        "tensile_strength": 505,
        "elongation": 40,
        "hardness_hv": 180,
        "youngs_modulus": 193,
        "density": 8.00,
        "thermal_conductivity": 16.2,
        "specific_heat": 500,
        "melting_point": 1400,
        "coefficient_thermal_expansion": 17.2,
        "electrical_resistivity": 72.0,
        "poissons_ratio": 0.29,
        "carbon_content": 0.08
    }
}

def main():
    st.markdown('<h1 class="main-header">⚙ Metallurgical Properties Calculator</h1>', 
                unsafe_allow_html=True)
    
    # Sidebar navigation
    st.sidebar.title("Navigation")
    calculation_type = st.sidebar.selectbox(
        "Select Calculation Type",
        [
            "Material Properties Database",
            "Mechanical Properties Calculator",
            "Phase Transformations",
            "Heat Treatment",
            "Corrosion Properties",
            "Material Selection",
            "Microstructure Analysis"
        ]
    )
    
    # Main content based on selection
    if calculation_type == "Material Properties Database":
        material_database_viewer()
    elif calculation_type == "Mechanical Properties Calculator":
        mechanical_properties_calculator()
    elif calculation_type == "Phase Transformations":
        phase_transformation_calculator()
    elif calculation_type == "Heat Treatment":
        heat_treatment_calculator()
    elif calculation_type == "Corrosion Properties":
        corrosion_properties_calculator()
    elif calculation_type == "Material Selection":
        material_selection_tool()
    elif calculation_type == "Microstructure Analysis":
        microstructure_analysis()

def material_database_viewer():
    st.header("📊 Material Properties Database")
    
    # Material selection
    selected_material = st.selectbox("Select Material", list(MATERIAL_DATABASE.keys()))
    
    if selected_material:
        material_data = MATERIAL_DATABASE[selected_material]
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader(f"Basic Properties - {selected_material}")
            st.markdown('<div class="result-box">', unsafe_allow_html=True)
            st.write(f"*Material Type:* {material_data['type']}")
            st.write(f"*Density:* {material_data['density']} g/cm³")
            st.write(f"*Melting Point:* {material_data['melting_point']}°C")
            st.write(f"*Young's Modulus:* {material_data['youngs_modulus']} GPa")
            st.write(f"*Poisson's Ratio:* {material_data['poissons_ratio']}")
            st.markdown('</div>', unsafe_allow_html=True)
            
            st.subheader("Mechanical Properties")
            st.markdown('<div class="result-box">', unsafe_allow_html=True)
            st.write(f"*Yield Strength:* {material_data['yield_strength']} MPa")
            st.write(f"*Tensile Strength:* {material_data['tensile_strength']} MPa")
            st.write(f"*Elongation:* {material_data['elongation']}%")
            st.write(f"*Hardness (HV):* {material_data['hardness_hv']}")
            st.markdown('</div>', unsafe_allow_html=True)
        
        with col2:
            st.subheader("Thermal Properties")
            st.markdown('<div class="result-box">', unsafe_allow_html=True)
            st.write(f"*Thermal Conductivity:* {material_data['thermal_conductivity']} W/m·K")
            st.write(f"*Specific Heat:* {material_data['specific_heat']} J/kg·K")
            st.write(f"*Coefficient of Thermal Expansion:* {material_data['coefficient_thermal_expansion']} ×10⁻⁶/°C")
            st.markdown('</div>', unsafe_allow_html=True)
            
            st.subheader("Electrical Properties")
            st.markdown('<div class="result-box">', unsafe_allow_html=True)
            st.write(f"*Electrical Resistivity:* {material_data['electrical_resistivity']} nΩ·m")
            st.markdown('</div>', unsafe_allow_html=True)
            
            # Calculate derived properties
            strength_ratio = material_data['yield_strength'] / material_data['tensile_strength']
            toughness_approx = (material_data['yield_strength'] + material_data['tensile_strength']) * material_data['elongation'] / 200
            specific_strength = material_data['tensile_strength'] / material_data['density']
            stiffness_index = material_data['youngs_modulus'] / material_data['density']
            
            st.subheader("Derived Properties")
            st.markdown('<div class="result-box">', unsafe_allow_html=True)
            st.write(f"*Strength Ratio (Y/T):* {strength_ratio:.3f}")
            st.write(f"*Approximate Toughness:* {toughness_approx:.1f} MPa·%")
            st.write(f"*Specific Strength:* {specific_strength:.1f} MPa·cm³/g")
            st.write(f"*Stiffness Index:* {stiffness_index:.1f} GPa·cm³/g")
            st.markdown('</div>', unsafe_allow_html=True)
        
        # Comparative charts
        st.subheader("Material Comparison")
        compare_materials = st.multiselect("Select materials to compare", 
                                         list(MATERIAL_DATABASE.keys()), 
                                         default=[selected_material])
        
        if compare_materials:
            properties_to_compare = st.selectbox("Select property to compare", 
                                              ["tensile_strength", "yield_strength", "elongation", 
                                               "hardness_hv", "density", "youngs_modulus"])
            
            materials = []
            values = []
            for material in compare_materials:
                materials.append(material)
                values.append(MATERIAL_DATABASE[material][properties_to_compare])
            
            fig = px.bar(x=materials, y=values, 
                        title=f"Comparison of {properties_to_compare.replace('_', ' ').title()}",
                        labels={"x": "Material", "y": properties_to_compare.replace('_', ' ').title()})
            st.plotly_chart(fig)

def mechanical_properties_calculator():
    st.header("Mechanical Properties Calculator")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Material Selection")
        selected_material = st.selectbox("Select Base Material", 
                                       list(MATERIAL_DATABASE.keys()),
                                       key="mech_mat_select")
        
        if selected_material:
            base_data = MATERIAL_DATABASE[selected_material]
            
            # Allow customization of base properties
            st.subheader("Customize Properties")
            yield_strength = st.number_input("Yield Strength (MPa)", 
                                           value=float(base_data['yield_strength']),
                                           min_value=0.0)
            tensile_strength = st.number_input("Tensile Strength (MPa)", 
                                            value=float(base_data['tensile_strength']),
                                            min_value=0.0)
            elongation = st.number_input("Elongation (%)", 
                                       value=float(base_data['elongation']),
                                       min_value=0.0)
            hardness = st.number_input("Hardness (HV)", 
                                     value=float(base_data['hardness_hv']),
                                     min_value=0.0)
            youngs_modulus = st.number_input("Young's Modulus (GPa)", 
                                           value=float(base_data['youngs_modulus']),
                                           min_value=0.0)
            density = st.number_input("Density (g/cm³)", 
                                    value=float(base_data['density']),
                                    min_value=0.0)
    
    with col2:
        st.subheader("Calculated Properties")
        
        # Calculate derived properties
        strength_ratio = yield_strength / tensile_strength if tensile_strength > 0 else 0
        toughness_approx = (yield_strength + tensile_strength) * elongation / 200
        specific_strength = tensile_strength / density if density > 0 else 0
        stiffness_index = youngs_modulus / density if density > 0 else 0
        
        st.markdown('<div class="result-box">', unsafe_allow_html=True)
        st.write(f"*Strength Ratio (Y/T):* {strength_ratio:.3f}")
        st.write(f"*Approximate Toughness:* {toughness_approx:.1f} MPa·%")
        st.write(f"*Specific Strength:* {specific_strength:.1f} MPa·cm³/g")
        st.write(f"*Stiffness Index:* {stiffness_index:.1f} GPa·cm³/g")
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Hardness conversions
        st.subheader("Hardness Conversions")
        hv_to_hb = hardness * 0.95  # Approximate conversion
        hv_to_hrc = 0.08 * hardness - 12.5  # Approximate conversion
        
        st.write(f"*Brinell Hardness (HB):* {hv_to_hb:.1f}")
        st.write(f"*Rockwell C (HRC):* {max(0, hv_to_hrc):.1f}")
        
        # Material classification
        st.subheader("Material Classification")
        if tensile_strength > 800:
            strength_class = "High Strength"
        elif tensile_strength > 400:
            strength_class = "Medium Strength"
        else:
            strength_class = "Low Strength"
        
        st.write(f"*Strength Class:* {strength_class}")
    
    # Stress-Strain Curve Simulation
    st.subheader("Stress-Strain Curve Simulation")
    
    # Generate synthetic stress-strain data
    strain_points = np.linspace(0, elongation/100 * 1.5, 100)
    stress_points = youngs_modulus * 1000 * strain_points  # Elastic region
    
    # Plastic region approximation
    plastic_start = yield_strength / (youngs_modulus * 1000)
    for i, strain in enumerate(strain_points):
        if strain > plastic_start:
            stress_points[i] = yield_strength + (tensile_strength - yield_strength) * \
                             ((strain - plastic_start) / (elongation/100 - plastic_start)) ** 0.5
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=strain_points, y=stress_points, mode='lines', name='Stress-Strain'))
    fig.add_trace(go.Scatter(x=[plastic_start], y=[yield_strength], 
                           mode='markers', name='Yield Point', marker=dict(size=10)))
    
    fig.update_layout(
        title="Engineering Stress-Strain Curve",
        xaxis_title="Strain",
        yaxis_title="Stress (MPa)",
        width=800
    )
    
    st.plotly_chart(fig)

def phase_transformation_calculator():
    st.header("Phase Transformation Calculator")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Steel Phase Transformation")
        
        # Material selection for phase transformations
        selected_material = st.selectbox("Select Steel Material", 
                                       [mat for mat in MATERIAL_DATABASE.keys() if "Steel" in mat],
                                       key="phase_mat_select")
        
        if selected_material:
            material_data = MATERIAL_DATABASE[selected_material]
            carbon_content = material_data.get('carbon_content', 0.2)
        else:
            carbon_content = 0.2
        
        carbon_content = st.slider("Carbon Content (%)", 0.02, 2.0, float(carbon_content), 0.01)
        alloy_type = st.selectbox("Alloy Type", ["Carbon Steel", "Low Alloy", "Stainless Steel"])
        
        # Calculate transformation temperatures (approximate)
        ae1 = 727 - 10.7 * carbon_content * 100  # Lower critical temperature
        ae3 = 910 - 203 * sqrt(carbon_content)   # Upper critical temperature
        ms = 539 - 423 * carbon_content - 30.4 * carbon_content**2  # Martensite start
        
        st.markdown('<div class="result-box">', unsafe_allow_html=True)
        st.write(f"*Ae1 Temperature:* {ae1:.1f}°C")
        st.write(f"*Ae3 Temperature:* {ae3:.1f}°C")
        st.write(f"*Ms Temperature:* {ms:.1f}°C")
        st.markdown('</div>', unsafe_allow_html=True)
        
        # TTT Diagram parameters
        st.subheader("TTT Diagram Parameters")
        incubation_time = st.number_input("Incubation Time (s)", min_value=0.1, value=10.0)
        pearlite_time = st.number_input("Pearlite Formation Time (s)", min_value=1.0, value=100.0)
        bainite_time = st.number_input("Bainite Formation Time (s)", min_value=1.0, value=50.0)
    
    with col2:
        st.subheader("Phase Fractions Calculator")
        
        cooling_rate = st.selectbox("Cooling Rate", ["Very Slow", "Slow", "Medium", "Fast", "Very Fast"])
        
        # Estimate phase fractions based on cooling rate
        cooling_rates = {"Very Slow": 0.1, "Slow": 1, "Medium": 10, "Fast": 100, "Very Fast": 1000}
        rate_factor = cooling_rates[cooling_rate]
        
        # Simplified phase fraction calculation
        ferrite_frac = max(0, 0.8 - 0.1 * rate_factor)
        pearlite_frac = max(0, 0.6 - 0.2 * rate_factor)
        bainite_frac = min(0.8, 0.3 + 0.1 * rate_factor)
        martensite_frac = min(1.0, 0.1 + 0.2 * rate_factor)
        
        # Normalize fractions
        total = ferrite_frac + pearlite_frac + bainite_frac + martensite_frac
        if total > 0:
            ferrite_frac /= total
            pearlite_frac /= total
            bainite_frac /= total
            martensite_frac /= total
        
        phases = ['Ferrite', 'Pearlite', 'Bainite', 'Martensite']
        fractions = [ferrite_frac, pearlite_frac, bainite_frac, martensite_frac]
        
        fig_pie = px.pie(values=fractions, names=phases, title="Estimated Phase Fractions")
        st.plotly_chart(fig_pie)
        
        # TTT Diagram simulation
        st.subheader("TTT Diagram Simulation")
        temperatures = np.linspace(200, 800, 50)
        times = incubation_time * np.exp(0.02 * (temperatures - 727))
        
        fig_ttt = go.Figure()
        fig_ttt.add_trace(go.Scatter(x=times, y=temperatures, mode='lines', name='Transformation Start'))
        fig_ttt.update_layout(title="Simplified TTT Diagram",
                            xaxis_title="Time (s)", yaxis_title="Temperature (°C)",
                            xaxis_type="log")
        st.plotly_chart(fig_ttt)

def heat_treatment_calculator():
    st.header("Heat Treatment Calculator")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.subheader("Material Selection")
        selected_material = st.selectbox("Select Material for Heat Treatment", 
                                       list(MATERIAL_DATABASE.keys()),
                                       key="heat_mat_select")
        
        if selected_material:
            material_data = MATERIAL_DATABASE[selected_material]
            initial_hardness = material_data['hardness_hv']
        else:
            initial_hardness = 200.0
            
        st.subheader("Austenitizing Parameters")
        aust_temp = st.number_input("Austenitizing Temperature (°C)", value=850.0)
        aust_time = st.number_input("Austenitizing Time (min)", value=30.0)
        
        st.subheader("Quenching Parameters")
        quench_medium = st.selectbox("Quenching Medium", 
                                   ["Water", "Oil", "Air", "Gas"])
        quench_rate = st.number_input("Quench Rate (°C/s)", value=50.0)
    
    with col2:
        st.subheader("Tempering Parameters")
        temper_temp = st.number_input("Tempering Temperature (°C)", value=200.0)
        temper_time = st.number_input("Tempering Time (min)", value=60.0)
        
        # Calculate tempered hardness (approximate)
        tempered_hardness = initial_hardness * (1 - 0.0005 * max(0, temper_temp - 200))
        
        st.markdown('<div class="result-box">', unsafe_allow_html=True)
        st.write(f"*Initial Hardness:* {initial_hardness:.1f} HV")
        st.write(f"*Estimated Tempered Hardness:* {tempered_hardness:.1f} HV")
        st.write(f"*Hardness Reduction:* {initial_hardness - tempered_hardness:.1f} HV")
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Quench severity calculation
        quench_severity = {"Water": 1.0, "Oil": 0.3, "Air": 0.02, "Gas": 0.1}
        severity = quench_severity[quench_medium]
        st.write(f"*Quench Severity (H-value):* {severity}")
    
    with col3:
        st.subheader("Jominy Test Calculator")
        jominy_distance = st.slider("Jominy Distance (mm)", 1, 50, 10)
        
        # Simplified Jominy curve calculation
        hardness_jominy = initial_hardness * np.exp(-0.05 * jominy_distance)
        
        st.write(f"*Hardness at {jominy_distance}mm:* {hardness_jominy:.1f} HV")
        
        # Generate Jominy curve
        distances = np.linspace(1, 50, 50)
        hardness_curve = initial_hardness * np.exp(-0.05 * distances)
        
        fig_jominy = go.Figure()
        fig_jominy.add_trace(go.Scatter(x=distances, y=hardness_curve, mode='lines'))
        fig_jominy.update_layout(title="Jominy Hardenability Curve",
                               xaxis_title="Distance from Quenched End (mm)",
                               yaxis_title="Hardness (HV)")
        st.plotly_chart(fig_jominy)

def corrosion_properties_calculator():
    st.header("Corrosion Properties Calculator")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Material Selection")
        selected_material = st.selectbox("Select Material for Corrosion Analysis", 
                                       list(MATERIAL_DATABASE.keys()),
                                       key="corr_mat_select")
        
        st.subheader("Corrosion Rate Calculation")
        
        weight_loss = st.number_input("Weight Loss (mg)", value=10.0)
        area = st.number_input("Exposed Area (cm²)", value=10.0)
        time = st.number_input("Exposure Time (days)", value=30.0)
        
        if selected_material:
            density = MATERIAL_DATABASE[selected_material]['density']
        else:
            density = 7.85
            
        # Corrosion rate calculation (mm/year)
        corrosion_rate = (weight_loss * 0.001) / (area * density * time) * 365 * 1000
        
        st.markdown('<div class="result-box">', unsafe_allow_html=True)
        st.write(f"*Corrosion Rate:* {corrosion_rate:.3f} mm/year")
        
        # Corrosion classification
        if corrosion_rate < 0.1:
            corr_class = "Excellent"
            st.success(f"*Corrosion Resistance:* {corr_class}")
        elif corrosion_rate < 0.5:
            corr_class = "Good"
            st.info(f"*Corrosion Resistance:* {corr_class}")
        elif corrosion_rate < 1.0:
            corr_class = "Fair"
            st.warning(f"*Corrosion Resistance:* {corr_class}")
        else:
            corr_class = "Poor"
            st.error(f"*Corrosion Resistance:* {corr_class}")
        
        st.markdown('</div>', unsafe_allow_html=True)
    
    with col2:
        st.subheader("Pitting Resistance")
        
        # For stainless steels, calculate PRE
        if selected_material and "Stainless" in selected_material:
            st.info("Stainless steel selected - calculating Pitting Resistance Equivalent")
            # Approximate PRE for common stainless steels
            pre_values = {
                "AISI 304 (Stainless Steel)": 25.0,
                "AISI 4140 (Alloy Steel)": 15.0  # Not stainless, but for example
            }
            pr_en = pre_values.get(selected_material, 20.0)
        else:
            pr_en = st.number_input("Pitting Resistance Equivalent (PRE)", value=20.0)
        
        if pr_en > 35:
            pitting_resistance = "Excellent"
            st.success(f"*Pitting Resistance:* {pitting_resistance}")
        elif pr_en > 25:
            pitting_resistance = "Good"
            st.info(f"*Pitting Resistance:* {pitting_resistance}")
        else:
            pitting_resistance = "Poor"
            st.warning(f"*Pitting Resistance:* {pitting_resistance}")
        
        # Environmental factors
        st.subheader("Environmental Factors")
        temperature = st.number_input("Environment Temperature (°C)", value=25.0)
        ph = st.slider("pH Level", 0.0, 14.0, 7.0)
        chloride_content = st.number_input("Chloride Content (ppm)", value=100.0)
        
        # Simple corrosion risk assessment
        risk_factor = (temperature / 50) + (abs(ph - 7) / 3.5) + (chloride_content / 1000)
        if risk_factor < 1:
            risk_level = "Low"
        elif risk_factor < 2:
            risk_level = "Medium"
        else:
            risk_level = "High"
        
        st.write(f"*Corrosion Risk Level:* {risk_level}")

def material_selection_tool():
    st.header("Material Selection Tool")
    
    # Extended material database for selection
    selection_db = {
        **MATERIAL_DATABASE,
        "Aluminum 6061-T6": {
            "type": "Aluminum Alloy", "YS": 276, "TS": 310, "Elong": 12, "Cost": 2.5, "density": 2.7
        },
        "Aluminum 7075-T6": {
            "type": "Aluminum Alloy", "YS": 503, "TS": 572, "Elong": 11, "Cost": 3.0, "density": 2.81
        },
        "Titanium Ti-6Al-4V": {
            "type": "Titanium Alloy", "YS": 830, "TS": 900, "Elong": 10, "Cost": 15.0, "density": 4.43
        }
    }
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Selection Criteria")
        min_yield = st.number_input("Minimum Yield Strength (MPa)", value=200.0)
        min_tensile = st.number_input("Minimum Tensile Strength (MPa)", value=300.0)
        min_elongation = st.number_input("Minimum Elongation (%)", value=10.0)
        max_cost = st.number_input("Maximum Cost Factor", value=10.0)
        max_density = st.number_input("Maximum Density (g/cm³)", value=10.0)
        
        material_type = st.multiselect("Material Types", 
                                     ["Carbon Steel", "Stainless Steel", "Aluminum Alloy", "Titanium Alloy"],
                                     default=["Carbon Steel", "Stainless Steel"])
    
    with col2:
        st.subheader("Selected Materials")
        
        suitable_materials = []
        for material, props in selection_db.items():
            if (props.get("yield_strength", props.get("YS", 0)) >= min_yield and 
                props.get("tensile_strength", props.get("TS", 0)) >= min_tensile and 
                props.get("elongation", props.get("Elong", 0)) >= min_elongation and 
                props.get("Cost", 5.0) <= max_cost and
                props.get("density", 8.0) <= max_density and
                props.get("type", "") in material_type):
                suitable_materials.append(material)
        
        if suitable_materials:
            st.success(f"Found {len(suitable_materials)} suitable materials:")
            for material in suitable_materials:
                props = selection_db[material]
                st.write(f"{material}")
                st.write(f"YS: {props.get('yield_strength', props.get('YS', 'N/A'))}MPa, "
                        f"TS: {props.get('tensile_strength', props.get('TS', 'N/A'))}MPa, "
                        f"Cost: {props.get('Cost', 'N/A')}x")
        else:
            st.warning("No materials meet the specified criteria")
        
        # Material comparison chart
        if suitable_materials:
            materials_list = suitable_materials[:5]  # Limit to 5 for clarity
            yield_strengths = [selection_db[m].get('yield_strength', selection_db[m].get('YS', 0)) for m in materials_list]
            tensile_strengths = [selection_db[m].get('tensile_strength', selection_db[m].get('TS', 0)) for m in materials_list]
            
            fig = go.Figure(data=[
                go.Bar(name='Yield Strength', x=materials_list, y=yield_strengths),
                go.Bar(name='Tensile Strength', x=materials_list, y=tensile_strengths)
            ])
            fig.update_layout(title="Selected Materials Strength Comparison")
            st.plotly_chart(fig)

def microstructure_analysis():
    st.header("Microstructure Analysis")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Grain Size Analysis")
        
        grain_size = st.number_input("Average Grain Size (μm)", value=50.0)
        magnification = st.number_input("Microscope Magnification", value=100.0)
        
        # ASTM grain size number
        astm_number = -6.64 * np.log10(grain_size * 0.001) - 3.29
        
        st.markdown('<div class="result-box">', unsafe_allow_html=True)
        st.write(f"*ASTM Grain Size Number:* {astm_number:.1f}")
        
        # Grain size classification
        if astm_number < 5:
            grain_class = "Coarse"
        elif astm_number < 8:
            grain_class = "Medium"
        else:
            grain_class = "Fine"
        
        st.write(f"*Grain Size Classification:* {grain_class}")
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Grain boundary analysis
        st.subheader("Grain Boundary Analysis")
        grain_boundary_length = st.number_input("Grain Boundary Length per unit area (μm/μm²)", value=0.5)
        grain_size_from_boundary = 2 / grain_boundary_length if grain_boundary_length > 0 else 0
        st.write(f"*Calculated Grain Size:* {grain_size_from_boundary:.1f} μm")
    
    with col2:
        st.subheader("Phase Analysis")
        phase1_frac = st.slider("Phase 1 Fraction (%)", 0, 100, 70)
        phase2_frac = 100 - phase1_frac
        
        phases = ['Phase 1', 'Phase 2']
        fractions = [phase1_frac, phase2_frac]
        
        fig = px.pie(values=fractions, names=phases, title="Phase Fraction Distribution")
        st.plotly_chart(fig)
        
        # Inclusion analysis
        st.subheader("Inclusion Analysis")
        inclusion_size = st.number_input("Maximum Inclusion Size (μm)", value=20.0)
        inclusion_density = st.number_input("Inclusion Density (per mm²)", value=100.0)
        
        # Cleanliness rating
        if inclusion_size < 10 and inclusion_density < 50:
            cleanliness = "Very Clean"
        elif inclusion_size < 20 and inclusion_density < 100:
            cleanliness = "Clean"
        else:
            cleanliness = "Contaminated"
        
        st.write(f"*Steel Cleanliness:* {cleanliness}")

# Add information section
def add_info_section():
    st.sidebar.markdown("---")
    st.sidebar.header("About")
    st.sidebar.info("""
    This application calculates various metallurgical properties including:
    - Material properties database with real data
    - Mechanical properties calculations
    - Phase transformations
    - Heat treatment parameters
    - Corrosion rates
    - Material selection
    - Microstructure analysis
    """)

if __name__ == "__main__":
    add_info_section()
    main()
