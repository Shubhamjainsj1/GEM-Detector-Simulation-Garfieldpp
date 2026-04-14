# 🚀 GEM Detector Simulation using Garfield++

## 📌 Overview
This project focuses on the simulation and analysis of **Gas Electron Multiplier (GEM) detectors** using **Garfield++**. The work includes studying electron transport, avalanche multiplication, gain behavior, and detector performance under different configurations.

---

## 🎯 Objectives
- Implement and validate **Single GEM simulation**
- Study **electron avalanche and gain behavior**
- Analyze **ion backflow**
- Investigate **effect of gas mixtures (Ar/CO₂)**
- Extend results towards **Triple GEM approximation**

---

## 🧠 Physics Background

Muon → Ionization → Primary electrons → GEM holes → Avalanche → Signal

- **Ionization**: Charged particles create electron-ion pairs  
- **Avalanche**: Strong electric field multiplies electrons  
- **Gain**: Ratio of output electrons to primary electrons  

---

## ⚙️ Technologies Used

- C++
- Garfield++
- ROOT
- Magboltz
- ANSYS field maps

---

## 🏗️ Project Structure

GEM-Detector-Simulation-Garfieldpp/
│── gem.C
│── gas_study.txt
│── ELIST.lis
│── NLIST.lis
│── MPLIST.lis
│── PRNSOL.lis
│── README.md

---

## 🧪 Simulation Details

### Gas Mixtures
- 90/10 (Ar/CO₂)
- 80/20
- 70/30
- 60/40
- 85/15

### Particle Setup
- Particle: Muon (μ⁻)
- Energy: 100 GeV
- Tool: TrackHeed

---

## 📊 Key Features

### ✅ Single GEM Simulation
- Electric field visualization  
- Electron drift and avalanche tracking  

### ✅ Gain Calculation

Gain = Total electrons after avalanche / Primary electrons

### ✅ Ion Backflow

IBF = Backflow ions / Total ions

### ✅ Gas Mixture Study
- Vector-based parameter scan  
- Comparison of gain and ion backflow  

### ✅ Triple GEM Approximation

G_total = (G1 + 12) × 28 × 28

---

## 📈 Results

### Observations
- Increasing CO₂:
  - Gain ↓
  - Ion backflow ↓
  - Stability ↑

- Increasing Argon:
  - Gain ↑
  - Ion backflow ↑

---

## 🛠️ How to Run

```bash
cd ~/garfieldpp/Examples/Gem
mkdir build
cd build
cmake ..
make
./Gem
