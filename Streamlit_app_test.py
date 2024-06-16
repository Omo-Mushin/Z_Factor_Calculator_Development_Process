
import logging



import streamlit as st
from functions import compute_Z_Hall_Yarborough, compute_Z_with_DPR, compute_Z_with_DAK

def main():
    st.title("Z-Factor Calculation")

    st.sidebar.header("Input Parameters")

    sp_gravity = st.sidebar.number_input("Enter specific gravity", min_value=0.0, format="%.4f")
    temperature = st.sidebar.number_input("Enter Temperature (Rankine)", min_value=0.0, format="%.2f")
    pressure = st.sidebar.number_input("Enter the pressure (psia)", min_value=0.0, format="%.2f")

    fluid_type = st.sidebar.selectbox("Select fluid type", ["NATURAL GAS", "CONDENSATE"])

    def Calculate_Z(sp_gravity, temperature, pressure, fluid_type):
        Z_HY = compute_Z_Hall_Yarborough(sp_gravity, temperature, pressure, fluid_type)
        Z_DAK = compute_Z_with_DAK(sp_gravity, temperature, pressure, fluid_type)
        # Z_DPR = compute_Z_with_DPR(sp_gravity, temperature, pressure, fluid_type)
        return Z_HY, Z_DAK

    if st.sidebar.button("Calculate"):
        try:
            Z_HY, Z_DAK, = Calculate_Z(sp_gravity, temperature, pressure, fluid_type)
            st.write(f'The Z-Factor using:')
            st.write(f'Hall-Yarborough : {Z_HY}')
            st.write(f'Dranchuk-Abu-Kaseem: {Z_DAK}')
            # st.write(f'Dranchuk-Purvis-Robinson: {Z_DPR}')
        except Exception as e:
            st.error(f'The error is {e}')


if __name__ == "__main__":
    main()
