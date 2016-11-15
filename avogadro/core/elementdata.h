#ifndef AVOGADRO_CORE_ELEMENTS_DATA
#define AVOGADRO_CORE_ELEMENTS_DATA

namespace Avogadro {
namespace Core {

unsigned char element_count = 119;

const char* element_symbols[] = {
  "Xx", "H", "He", "Li", "Be", "B", "C", "N", "O", "F",
  "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K",
  "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
  "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
  "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
  "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
  "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
  "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
  "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
  "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
  "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };

const char* element_names[] = {
  "Dummy", "Hydrogen", "Helium", "Lithium", "Beryllium",
  "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine",
  "Neon", "Sodium", "Magnesium", "Aluminium", "Silicon",
  "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium",
  "Calcium", "Scandium", "Titanium", "Vanadium", "Chromium",
  "Manganese", "Iron", "Cobalt", "Nickel", "Copper",
  "Zinc", "Gallium", "Germanium", "Arsenic", "Selenium",
  "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium",
  "Zirconium", "Niobium", "Molybdenum", "Technetium", "Ruthenium",
  "Rhodium", "Palladium", "Silver", "Cadmium", "Indium",
  "Tin", "Antimony", "Tellurium", "Iodine", "Xenon",
  "Caesium", "Barium", "Lanthanum", "Cerium", "Praseodymium",
  "Neodymium", "Promethium", "Samarium", "Europium", "Gadolinium",
  "Terbium", "Dysprosium", "Holmium", "Erbium", "Thulium",
  "Ytterbium", "Lutetium", "Hafnium", "Tantalum", "Tungsten",
  "Rhenium", "Osmium", "Iridium", "Platinum", "Gold",
  "Mercury", "Thallium", "Lead", "Bismuth", "Polonium",
  "Astatine", "Radon", "Francium", "Radium", "Actinium",
  "Thorium", "Protactinium", "Uranium", "Neptunium", "Plutonium",
  "Americium", "Curium", "Berkelium", "Californium", "Einsteinium",
  "Fermium", "Mendelevium", "Nobelium", "Lawrencium", "Rutherfordium",
  "Dubnium", "Seaborgium", "Bohrium", "Hassium", "Meitnerium",
  "Darmstadtium", "Roentgenium", "Copernicium", "Nihonium", "Flerovium",
  "Moscovium", "Livermorium", "Tennessine", "Oganesson" };

double element_masses[] = {
  // from IUPAC http://www.chem.qmul.ac.uk/iupac/AtWt/
  // (2015 set, updated from 2013)
  0, 1.00784, 4.0026,
  6.938, 9.01218, 10.806, 12.011, 14.006, 15.9994, 18.9984, 20.1797,
  22.9898, 24.305, 26.9815,
  28.0855, 30.9738, 32.065, 35.453, 39.948, 39.0983, 40.078,
  44.9559, 47.867, 50.9415, 51.9961, 54.938, 55.845, 58.9332,
  58.6934, 63.546, 65.38, 69.723, 72.64, 74.9216, 78.971,
  79.904, 83.798, 85.4678, 87.62, 88.9058, 91.224, 92.9064,
  95.95, 97, 101.07, 102.9055, 106.42, 107.8682, 112.414,
  114.818, 118.71, 121.76, 127.6, 126.9045, 131.293, 132.9055,
  137.327, 138.9055, 140.116, 140.9077, 144.242, 145, 150.36,
  151.964, 157.25, 158.9253, 162.5, 164.9303, 167.259, 168.9342,
  173.045, 174.9668, 178.49, 180.9479, 183.84, 186.207, 190.23,
  192.217, 195.084, 196.9666, 200.592, 204.38, 207.2, 208.9804,
  209, 210, 222, 223, 226, 227, 232.0377,
  231.0358, 238.0289, 237, 244, 243, 247, 247,
  251, 252, 257, 258, 259, 262, 267,
  270, 269, 270, 270, 278, 281, 281,
  285, 286, 289, 289, 293, 293, 294 };

double element_VDW[] = {
  // From Alvarez doi: 10.1039/C3DT50599E
  // Dalton Trans., 2013,42, 8617-8636
  // Dummy, 1st row
  0.69, 1.2, 1.43,
  // 2nd row (Li..Ne)
  2.12, 1.98, 1.91, 1.77, 1.66, 1.50, 1.46, 1.58,
  // 3rd row (Na .. Ar)
  2.50, 2.51, 2.25, 2.19, 1.90, 1.89, 1.82, 1.83,
  // 4th row (K, Ca)
  2.73, 2.62,
  // 1st row TM (Sc.. Zn)
  2.58, 2.46, 2.42, 2.45, 2.45, 2.44, 2.40, 2.40, 2.38, 2.39,
  // 4th row p-block (Ga .. Kr)
  2.32, 2.29, 1.88, 1.82, 1.86, 2.25,
  // 5th row Rb, Sr
  3.21, 2.84,
  // 2nd row TM (Y .. Cd)
  2.75, 2.52, 2.56, 2.45, 2.44, 2.46, 2.44, 2.15, 2.53, 2.49,
  // 5th row p-block (Sn .. Xe)
  2.43, 2.42, 2.47, 1.99, 2.04, 2.06,
  // 6th row Cs, Ba
  3.48, 3.03,
  // Lanthanides (La..Gd)
  2.98, 2.88, 2.92, 2.95, 2.90, 2.87, 2.83,
  // Lanthanides (Tb..Yb)
  2.79, 2.87, 2.81, 2.83, 2.79, 2.80,
  // 3rd row TM (Lu..Hg)
  2.74, 2.63, 2.53, 2.57, 2.49, 2.48, 2.41, 2.29, 2.32, 2.45,
  // 6th row p-block (Tl.. Bi)
  // 2.5 is a default here
  2.47, 2.60, 2.54, 2.5, 2.5, 2.5,
  // 7th row
  // 2.5 is a default here
  2.5, 2.5,
  // Actinides
  2.8, 2.93, 2.88, 2.71, 2.82, 2.81, 2.83,
  3.05, 3.38, 3.05, 3., 3., 3., 3.,
  // Trans-actinides
  3., 3., 3., 3., 3., 3., 3., 3., 3., 3.,
  // 7th row p-block
  3., 3., 3., 3., 3., 3.,
};


double element_covalent[] = {
  // From Pyykko doi: 10.1002/chem.200800987
  // Dummy, 1st row
  0.18, 0.32, 0.46,
  // 2nd row
  1.33, 1.02, 0.85, 0.75, 0.71, 0.63, 0.64, 0.67,
  // 3rd row
  1.55, 1.39, 1.26, 1.16, 1.11, 1.03, 0.99, 0.96,
  // 4th row K, Ca
  1.96,  1.71,
  // 1st row TM (Sc.. Zn)
  1.48, 1.36, 1.34, 1.22, 1.19, 1.16, 1.11, 1.10, 1.12, 1.18,
  // 4th row p-block (Ga..Kr)
 1.24, 1.21, 1.21, 1.16, 1.14, 1.17,
  // 5th row Rb, Sr
 2.10, 1.85,
 // 2nd row TM (Y..Cd)
  1.63, 1.54, 1.47, 1.38, 1.28, 1.25, 1.25, 1.20, 1.28, 1.36,
  // 5th row p-block (In..Xe)
  1.42, 1.40, 1.40, 1.36, 1.33, 1.31,
  // 6th row Cs, Ba
  2.32, 1.96,
  // Lanthanides La..Gd
  1.80, 1.63, 1.76, 1.74, 1.73, 1.72, 1.68,
  // Lanthanides Tb..Yb
  1.69, 1.68, 1.67, 1.66, 1.65, 1.64, 1.70,
  // 3rd row TM (Lu..Hg)
  1.62, 1.52, 1.46, 1.37, 1.31, 1.29, 1.22, 1.23, 1.24, 1.33,
  // 6th row p-block (Tl..Rn)
  1.44, 1.44, 1.51, 1.45, 1.47, 1.42,
  // 7th row Fr, Ra
  2.23, 2.01,
  // Actinides (Ac.. Am)
  1.86, 1.75, 1.69, 1.70, 1.71, 1.72, 1.66,
  // Actinides (Cm..No)
  1.66, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76,
  // Trans-actinides
  1.61, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28, 1.21, 1.22,
  1.36, 1.43, 1.62, 1.75, 1.65, 1.57
};

unsigned char element_color[][3] = {
  {17, 127, 178}, {255, 255, 255}, {216, 255, 255},
  {204, 127, 255}, {193, 255, 0}, {255, 181, 181},
  {127, 127, 127}, {12, 12, 255}, {255, 12, 12},
  {178, 255, 255}, {178, 226, 244}, {170, 91, 242},
  {137, 255, 0}, {191, 165, 165}, {127, 153, 153},
  {255, 127, 0}, {255, 255, 48}, {30, 239, 30},
  {127, 209, 226}, {142, 63, 211}, {61, 255, 0},
  {229, 229, 229}, {191, 193, 198}, {165, 165, 170},
  {137, 153, 198}, {155, 122, 198}, {127, 122, 198},
  {112, 122, 198}, {91, 122, 193}, {255, 122, 96},
  {124, 127, 175}, {193, 142, 142}, {102, 142, 142},
  {188, 127, 226}, {255, 160, 0}, {165, 40, 40},
  {91, 183, 209}, {112, 45, 175}, {0, 255, 0},
  {147, 255, 255}, {147, 224, 224}, {114, 193, 201},
  {84, 181, 181}, {58, 158, 158}, {35, 142, 142},
  {10, 124, 140}, {0, 104, 132}, {224, 224, 255},
  {255, 216, 142}, {165, 117, 114}, {102, 127, 127},
  {158, 99, 181}, {211, 122, 0}, {147, 0, 147},
  {66, 158, 175}, {86, 22, 142}, {0, 201, 0},
  {112, 211, 255}, {255, 255, 198}, {216, 255, 198},
  {198, 255, 198}, {163, 255, 198}, {142, 255, 198},
  {96, 255, 198}, {68, 255, 198}, {48, 255, 198},
  {30, 255, 198}, {0, 255, 155}, {0, 229, 117},
  {0, 211, 81}, {0, 191, 56}, {0, 170, 35},
  {76, 193, 255}, {76, 165, 255}, {33, 147, 214},
  {38, 124, 170}, {38, 102, 150}, {22, 84, 135},
  {244, 237, 209}, {204, 209, 30}, {181, 181, 193},
  {165, 84, 76}, {86, 89, 96}, {158, 79, 181},
  {170, 91, 0}, {117, 79, 68}, {66, 130, 150},
  {66, 0, 102}, {0, 124, 0}, {112, 170, 249},
  {0, 186, 255}, {0, 160, 255}, {0, 142, 255},
  {0, 127, 255}, {0, 107, 255}, {84, 91, 242},
  {119, 91, 226}, {137, 79, 226}, {160, 53, 211},
  {178, 30, 211}, {178, 30, 186}, {178, 12, 165},
  {188, 12, 135}, {198, 0, 102}, {204, 0, 89},
  {209, 0, 79}, {216, 0, 68}, {224, 0, 56},
  {229, 0, 45}, {232, 0, 38}, {234, 0, 35},
  {237, 0, 33}, {239, 0, 30}, {242, 0, 28},
  {244, 0, 25}, {247, 0, 22}, {249, 0, 20},
  {252, 0, 17}, {255, 0, 15} };

}
}

#endif
