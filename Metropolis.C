/// \file
/// \ingroup tutorial_graphics
/// Macro illustrating how to animate a picture using a Timer.
///
/// \macro_code
///
/// \author Rene Brun

#include <algorithm>
#include <random>
#include <condition_variable>

#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH2I.h"
#include "TTimer.h"

// critical temperature for 2D Ising model, based on http://web.mit.edu/krish_s/www/files/ising_Model.pdf
const Double_t Tc = 2.269; // using k_B = 1
const UInt_t initial_size = 100; // length of (square) lattice - it will have periodic boundaries!
std::vector<std::vector<Int_t>> s; // 2D "array" of spins (= +1 or -1, each), vector -> can be resized
Double_t T = Tc; // temperature of Ising model
Double_t B = 0.; // magnetic field
Double_t J = 1.; // coupling strength

std::vector<UInt_t> forbidden_i;
std::vector<UInt_t> forbidden_j;
UInt_t max_non_adjacent = 0;
const UInt_t k_max_inital = 15;
UInt_t k_max = k_max_inital;

TH2I *f2;
TCanvas* c1;
TTimer* timer;

std::mt19937 rng(1);    // random-number engine used (Mersenne-Twister in this case)
std::uniform_real_distribution<Double_t> uni(0.,1.); // guaranteed unbiased
std::uniform_int_distribution<UInt_t> uni_int(0,99);

std::mutex mtx;
std::condition_variable cv;

// Below are variables used for calculating thermodynamic properties

std::vector<Double_t> m_vec;
std::vector<Double_t> E_vec;
std::vector<Double_t> C_vec;
std::vector<Double_t> m2_vec; // used for susceptibilit

std::vector<Double_t> m_val_vec; // magnetisation
std::vector<Double_t> E_val_vec; // Energy (expectation value)
std::vector<Double_t> C_val_vec; // Heat capacity
std::vector<Double_t> X_val_vec; // magnetic susceptibility

UInt_t n = 5000;
UInt_t counter = 0;

Bool_t kCalcT = kFALSE;
Bool_t kCalcB = kFALSE;

Double_t m_val;
Double_t E_val;
Double_t X_val;

const UInt_t NdT = 200;

const Double_t dt = 0.001;
const Double_t min_T = Tc - dt * ((Double_t) NdT) / 2.;
const Double_t max_T = Tc + dt * ((Double_t) NdT) / 2.;
const UInt_t N_steps_relax = 1000000;

const Double_t dB = 0.01;
const Double_t min_B = -2. * dB;
const Double_t max_B = -min_B;

//File for saving results (running in batch mode speeds up calc())
TFile* f;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void setB(Double_t val){
    B = val;
    return;
}

void setT(Double_t val){
    if(val > 0.) T = val;
    else std::cout << "Cannot set T <= 0." << std::endl;
    return;
}

void setJ(Double_t val){
    J = val;
    return;
}

void GetMComp(){
    
    UInt_t size = s.size();
    Double_t sum = 0;
    
    for (UInt_t i = 0; i < size; i++) {
        for (UInt_t j = 0; j < size; j++) {
            sum += (Double_t) s.at(i).at(j);
        }
    }
    
    m_vec.push_back(sum);
    return;
}

Double_t GetM(){
    
    Double_t size = (Double_t) s.size();
    Double_t N = (size * size);
    
    Double_t sum = 0;
    for (UInt_t i = 0; i < m_vec.size(); i++) {
        sum += m_vec.at(i);
    }
    
    sum /= (Double_t) m_vec.size();
    
    sum /= N;
    return sum;
}

void ClearM(){
    m_vec.clear();
    return;
}

void GetEComp(){
    
    UInt_t size = s.size();
    Double_t E = 0;
    
    Double_t sum = 0.;
    
    if(B != 0.){
    for (UInt_t i = 0; i < size; i++) {
        for (UInt_t j = 0; j < size; j++) {
            sum += (Double_t) s.at(i).at(j);
        }
    }
    
    E += -B * sum;
    }
    
    
    // now sum over all pairs. We'll count each twice, so we'll need to divide by 2 at the end
    sum = 0.;
    
    Double_t n1, n2, n3, n4;
    
    for (UInt_t i = 0; i < size; i++) {
        for (UInt_t j = 0; j < size; j++) {
            
            if(j == 0) n1 = (Double_t) s.at(i).at(size-1);
            else n1 = (Double_t) s.at(i).at(j-1);
            
            if(i == size-1) n2 = (Double_t) s.at(0).at(j);
            else n2 = (Double_t) s.at(i+1).at(j);
            
            if(j == size-1) n3 = (Double_t) s.at(i).at(0);
            else n3 = (Double_t) s.at(i).at(j+1);
            
            if(i == 0) n4 = (Double_t) s.at(size-1).at(j);
            else n4 = (Double_t) s.at(i-1).at(j);
            
            sum += ((Double_t) s.at(i).at(j)) * (n1 + n2 + n3 + n4);
        }
    }
    sum /= 2.;
    
    E+= -J * sum;
    
    E_vec.push_back(E);
    return;
}

Double_t GetE(){
    
    Double_t sum = 0.;
    Double_t size = (Double_t) E_vec.size();
    for (UInt_t i = 0; i < E_vec.size(); i++) {
        sum += E_vec.at(i);
    }
    
    sum /= size;
    return sum;
}

void ClearE(){
    
    E_vec.clear();
    return;
}

Double_t GetX(){
    
    // Method: We have a bunch of m's and know dB. With length l, find l-1 "simple" derivatives, then average
    UInt_t length = m2_vec.size();
    
    Double_t val1, val2, sum = 0.;
    for (UInt_t i = 0; i < length - 1; i++) {
        val1 = m2_vec.at(i);
        val2 = m2_vec.at(i+1);
        sum += (val2 - val1) / dB;
    }
    
    return sum / ((Double_t)(length - 1));
}

//smooth out our m vector
void CleanM(){
    for (UInt_t i = 1; i < m_val_vec.size()-1; i++) m_val_vec.at(i) = (m_val_vec.at(i-1) + m_val_vec.at(i) + m_val_vec.at(i+1))/3.;
    return;
}

//smooth out our E vector
void CleanE(){
    for (UInt_t i = 1; i < E_val_vec.size()-1; i++) E_val_vec.at(i) = (E_val_vec.at(i-1) + E_val_vec.at(i) + E_val_vec.at(i+1))/3.;
    return;
}

void GetC(){
    
    // C = d<E>/dT. Therefore we need the E_val_vec, which contains <E> at different values of T.
    // This function will benefit from a choice of a small dt.
    
    UInt_t length = E_val_vec.size();
    if(length == 0) return;
    
    C_val_vec.push_back(0.); // first entry
    
    Double_t E_before, E, E_after, deriv1, deriv2, val;
    for (UInt_t i = 1; i < length - 1; i++) {
        
        E_before = E_val_vec.at(i-1);
        E = E_val_vec.at(i);
        E_after = E_val_vec.at(i+1);
        
        deriv1 = (E - E_before) / dt;
        deriv2 = (E_after - E) / dt;
        val = (deriv1 + deriv2) / 2.;
        C_val_vec.push_back(val);
    }
    C_val_vec.push_back(0.); // last entry -- now C_val_vec has same length as other value vectors (e.g. E_val_vec)
    
    //let's linearize the ends of the C_val_vec.
    C_val_vec.at(0) = C_val_vec.at(1) - (C_val_vec.at(2) - C_val_vec.at(1));
    C_val_vec.at(C_val_vec.size()-1) = C_val_vec.at(C_val_vec.size()-2) + (C_val_vec.at(C_val_vec.size()-2) - C_val_vec.at(C_val_vec.size()-3));
    return;
}

void init(UInt_t size = initial_size, Bool_t start = kTRUE){
    
    m_val_vec.clear();
    E_val_vec.clear();
    C_val_vec.clear();

    if(start){
        s.clear();
        std::vector<Int_t> s_row;
        
        for (UInt_t i = 0; i < size; i++) {
            
            Int_t val;
            
            for (UInt_t j = 0; j < size; j++) {
                
                if(uni(rng) > 0.5) val = -1;
                else val = 1;
                s_row.push_back(val);
            }
            s.push_back(s_row);
            s_row.clear();
        }
    }
    //~~~~~~~~~~~
    
    gStyle->SetOptStat(0);
    gStyle->SetCanvasPreferGL(true);
    //    gStyle->SetFrameFillColor(42);
    c1 = new TCanvas("c1");
    f2 = new TH2I("MMC","MMC", size,0,size,size,0,size);
    
    for(UInt_t i = 0; i < size; i++){
        for (UInt_t j = 0; j < size; j++) {
            f2->SetBinContent(i+1,j+1, s.at(i).at(j));
        }
    }
    f2->SetMaximum(1);
    f2->SetMinimum(-1);
    f2->Draw("COLZ");
    
    if(size % 2 == 0){
        
        max_non_adjacent = size / 2;
        max_non_adjacent *= size;
    }
    else{
        
        max_non_adjacent = (size + 1) / 2;
        max_non_adjacent = max_non_adjacent * max_non_adjacent + (max_non_adjacent - 1) * (max_non_adjacent - 1);
    }
    
    k_max = k_max_inital;
    if(max_non_adjacent < k_max_inital) k_max = max_non_adjacent;
    
    return;
}

// calculate change in energy given a spin flip
Double_t deltaE(UInt_t i, UInt_t j, UInt_t size){
    
    /*
     * values of 4 neighbors
     *  ^    -<-i->+
     *  |       n1
     *  j   n4      n2
     *  |       n3
     *  v
     */
    Int_t n1, n2, n3, n4;
    
    if(j == 0) n1 = s.at(i).at(size-1);
    else n1 = s.at(i).at(j-1);
    
    if(i == size-1) n2 = s.at(0).at(j);
    else n2 = s.at(i+1).at(j);
    
    if(j == size-1) n3 = s.at(i).at(0);
    else n3 = s.at(i).at(j+1);
    
    if(i == 0) n4 = s.at(size-1).at(j);
    else n4 = s.at(i-1).at(j);
    
    Int_t flipped = -s.at(i).at(j);
    
    Double_t dE = -B * ((Double_t) flipped) - J * ((Double_t) flipped) * ((Double_t) (n1 + n2 + n3 + n4));
    return dE;
}

// simulate one step
std::vector<UInt_t> step(UInt_t size){
    
    UInt_t arr [2];
    
    Bool_t cont = kFALSE;
    UInt_t i = 0;
    UInt_t j = 0;
    while(!cont){
        
        i = uni_int(rng) % size; // need the mod in case of rescaling
        j = uni_int(rng) % size;
        
        cont = kTRUE;
        
        // repeat if we pick a forbidden spin (i.e. one adjacent to one we're simultaneously changing)
        for(UInt_t k = 0; k < forbidden_i.size(); k++){
            if(i == forbidden_i.at(k) && j == forbidden_j.at(k)){
                cont = kFALSE;
                break;
            }
        }
    }
    
    Double_t dE = deltaE(i, j, size);
    
    if (dE < 0) s.at(i).at(j) *= -1;
    else if (uni(rng) < exp(-dE / T)) s.at(i).at(j) *= -1;
    
    std::vector<UInt_t> vec;
    vec.push_back(i);
    vec.push_back(j);
    
    return vec;
}

// rescale - use "majority rules"
void rg(UInt_t power = 1){
    
    for(UInt_t k = 0; k < power; k++){
        
        UInt_t size_old = s.size();
        if (size_old <= 2) return;
        UInt_t size_new = size_old / 2;
        if(size_new % 2 !=0) size_new --;
        
        std::vector<std::vector<Int_t>> s_new;
        
        Int_t val;
        std::vector<Int_t> s_new_row;
        for (UInt_t i = 0; i < size_old; i += 2) {
            for (UInt_t j = 0; j < size_old; j += 2) {
                
                val = s.at(i).at(j) + s.at(i).at(j+1) + s.at(i+1).at(j) + s.at(i+1).at(j+1);
                if(val > 0) val = 1;
                else if (val < 0) val = -1;
                else if (val == 0){
                    if(uni(rng)>0.5) val = 1;
                    else val = -1;
                }
                
                s_new_row.push_back(val);
            }
            s_new.push_back(s_new_row);
            s_new_row.clear();
        }
        
        s.clear();
        
        for (UInt_t i = 0; i < size_new; i++) {
            s.push_back(s_new.at(i));
        }
        
    }
    
    delete f2;
    delete c1;
    init(s.size(), kFALSE);
    return;
}

void run(UInt_t tStep = 10){
    timer = new TTimer(tStep);
    timer->SetCommand("Animate()");
    timer->TurnOn();
}

void stop(){
    timer->TurnOff();
}

void calc(Double_t temp = min_T, Double_t bfield = 0., Bool_t bool_temp = kTRUE){
    
    //cleanup
    ClearM();
    ClearE();
    counter = 0;
    
    setT(temp);
    setB(bfield);
    
    //vars for copy of Animate()
    UInt_t size;
    std::vector<UInt_t> vec;
    // fast forward animation after changing temp/B-field (simplified Animate())
    for (UInt_t h = 0; h < N_steps_relax; h++) {
        forbidden_i.clear();
        forbidden_j.clear();
        size = s.size();
        for (UInt_t k = 0; k < k_max; k++) {
            vec = step(size);
            UInt_t i = vec.at(0);
            UInt_t j = vec.at(1);
            f2->SetBinContent(i+1,j+1,s.at(i).at(j));
        }
    }
    
    if(bool_temp) kCalcT = kTRUE;
    else kCalcB = kTRUE;
    return;
}

void PlotM(UInt_t option = 0){
    
    if(m_val_vec.size() < NdT){
        std::cout << "Ending. NdT =      " << NdT << std::endl;
        std::cout << "m_val_vec.size() = " << m_val_vec.size() << std::endl;
        return;
    }
    
    Double_t tau = 0.;
    TCanvas* c2 = new TCanvas("c2", "c2", 600, 400);
    TH1F* hist;
    
    if(option == 0) hist = new TH1F("hM", "m(#tau);#tau", NdT, min_T - Tc, max_T - Tc);
    else hist = new TH1F("hM", "m(T);T", NdT, min_T, max_T);
    for (UInt_t i = 0; i < NdT; i++) {
        hist->SetBinContent(i, m_val_vec.at(i));
    }
    
    //our histogram may be noisy, let's smooth the data with a 3-point spline (I think that's what this is?)
    Double_t binb, bina;
    for(UInt_t i = 1; i < NdT-1; i++){
        binb = hist->GetBinContent(i-1);
        bina = hist->GetBinContent(i+1);
        hist->SetBinContent(i, (binb + bina + hist->GetBinContent(i))/3.);
    }
    
    c2->cd();
    hist->Draw("HIST");
    
    //warning: bad if file is not opened via sim()
    f->cd();
    hist->Write("hM");
    c2->Write("c2");
    return;
}

void PlotC(UInt_t option = 0){
    
    Double_t N_spins = (Double_t)(s.size() * s.size());
    if(C_val_vec.size() < NdT){
        std::cout << "Ending. NdT =      " << NdT << std::endl;
        std::cout << "C_val_vec.size() = " << C_val_vec.size() << std::endl;
        return;
    }
    
    Double_t tau = 0.;
    TCanvas* c3 = new TCanvas("c3", "c3", 600, 400);
    TH1F* hist;
    
    if(option == 0) hist = new TH1F("hC", "c(#tau);#tau", NdT, min_T - Tc, max_T - Tc);
    else hist = new TH1F("hC", "c(T);T", NdT, min_T, max_T);
    for (UInt_t i = 0; i < NdT; i++) {
        hist->SetBinContent(i, C_val_vec.at(i));
    }
    //our histogram may be noisy, let's smooth the data with a 3-point spline (I think that's what this is?)
    Double_t binb, bina;
    for(UInt_t i = 1; i < NdT-1; i++){
        binb = hist->GetBinContent(i-1);
        bina = hist->GetBinContent(i+1);
        hist->SetBinContent(i, (binb + bina + hist->GetBinContent(i))/3.);
    }
    
    //get histogram of c = C/N
    hist->Scale(1./N_spins);
    
    c3->cd();
    hist->Draw("HIST");
    
    //warning: bad if file is not opened via sim()
    f->cd();
    hist->Write("hC");
    c3->Write("c3");
    return;
}

void PlotX(UInt_t option = 0){
    
//    if(C_val_vec.size() < NdT){
//        std::cout << "Ending. NdT =      " << NdT << std::endl;
//        std::cout << "X_val_vec.size() = " << X_val_vec.size() << std::endl;
//        return;
//    }
    
    Double_t tau = 0.;
    TCanvas* c4 = new TCanvas("c4", "c4", 600, 400);
    TH1F* hist;
    
    if(option == 0) hist = new TH1F("hX", "#xi(#tau);#tau", NdT, min_T - Tc, max_T - Tc);
    else hist = new TH1F("hX", "#xi(T);T", NdT, min_T, max_T);
    for (UInt_t i = 0; i < NdT; i++) {
        hist->SetBinContent(i, X_val_vec.at(i));
    }
    //our histogram may be noisy, let's smooth the data with a 3-point spline (I think that's what this is?)
    Double_t binb, bina;
    for(UInt_t i = 1; i < NdT-1; i++){
        binb = hist->GetBinContent(i-1);
        bina = hist->GetBinContent(i+1);
        hist->SetBinContent(i, (binb + bina + hist->GetBinContent(i))/3.);
    }
    
    c4->cd();
    hist->Draw("HIST");
    
    //warning: bad if file is not opened via sim()
    f->cd();
    hist->Write("hX");
    c4->Write("c4");
    return;
}

void Animate(Bool_t update = kTRUE)
{
    
    forbidden_i.clear();
    forbidden_j.clear();
    UInt_t size = s.size();
    std::vector<UInt_t> vec;
    
    for (UInt_t k = 0; k < k_max; k++) {
        
        vec = step(size);
        //just in case the canvas has been deleted
        if (!gROOT->GetListOfCanvases()->FindObject("c1")) return;
        UInt_t i = vec.at(0);
        UInt_t j = vec.at(1);
        f2->SetBinContent(i+1,j+1,s.at(i).at(j));
        
        if(kCalcT || kCalcB){
            if(counter < n){
                GetMComp();
                GetEComp();
                counter++;
            }
            else{
                m_val = GetM();
                // calculating m, E
                if(kCalcT){
                    E_val = GetE();
                    m_val_vec.push_back(m_val);
                    E_val_vec.push_back(E_val);
                    if(T >= max_T){
                        kCalcT = kFALSE;
                        std::cout << "Calc'd - T." << std::endl;
                        PlotM(1);
                        GetC();
                        PlotC(1);
                    }
                    else calc(T + dt);
                }
                // calculating X
                else{
                    m2_vec.push_back(m_val);
                    if(B >= max_B){
                        X_val = GetX();
                        X_val_vec.push_back(X_val);
                        m2_vec.clear();
                        if(T >= max_T){
                            kCalcB = kFALSE;
                            std::cout << "Calc'd - B." << std::endl;
                            PlotX(1);
                        }
                        else calc(T + dt, min_B, kFALSE);
                    }
                    else calc(T, B + dB, kFALSE);
                }
            }
        }
        if(update){
            gPad->Modified();
            gPad->Update();
        }
    }
    return;
}

void sim(UInt_t option = 0, TString root_name = "calc.root"){
    
    f = new TFile(root_name, "UPDATE");
    init(100);
    run(1);
    if(option == 0) calc(min_T, 0., kTRUE);
    else calc(min_T, min_B, kFALSE);
    return;
}

// remember to close the file whenever plotting is done
void cf(){
    f->Close();
}

