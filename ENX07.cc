// ENX07_extended.cc

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserEventAction.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserRunAction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4GenericIon.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4Event.hh"

#include <filesystem>
#include <fstream>
#include <map>

// namespace ENX07 {
  
// 1) EventAction — акумуляція загального Edep
class EventAction : public G4UserEventAction {
public:
  void BeginOfEventAction(const G4Event* evt) override {
    edepTotal = 0.;
  }
  void EndOfEventAction(const G4Event*) override {
    if (edepTotal>0.)
      G4cout << "Event total Edep = " << edepTotal/keV << " keV\n";
  }
  void AddEdep(G4double ed) { edepTotal += ed; }
private:
  G4double edepTotal = 0.;
};

// 2) Sensitive detector — реєстрація енергії та процесів
class SpectrumDetector : public G4VSensitiveDetector {
public:
  SpectrumDetector(const G4String& name)
    : G4VSensitiveDetector(name)
  {
    std::filesystem::create_directories("logs");
    fout.open("logs/detected_energies.txt");
  }
  ~SpectrumDetector() override {
    if (fout.is_open()) fout.close();
  }

  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
    auto* track = step->GetTrack();
    if (track->GetDefinition() != G4Gamma::GammaDefinition())
      return false;

    G4double edep = step->GetTotalEnergyDeposit();
    G4double kinE = step->GetPreStepPoint()->GetKineticEnergy();
    if (edep<=0.) return false;

    // 1) Назва процесу, що згенерував крок
    auto* post = step->GetPostStepPoint();
    auto* proc = post->GetProcessDefinedStep();
    G4String name = proc
        ? proc->GetProcessName()
        : G4String("unknown");

    // 2) Ідентифікація ID для гістограми енергії
    static std::map<G4String,int> procID = {
      {"phot", 4},    // фотоефект
      {"compt",5},    // Комптон
      {"conv", 6}     // народження пари
    };
    int id = procID.count(name) ? procID[name] : 7;  // 7 — інші

    // 3) Заповнюємо гістограми
    auto* man = G4AnalysisManager::Instance();
    man->FillH1(1, edep/keV);   // DetectedGamma (усі)
    man->FillH1(2, edep/keV);   // EdepFine
    man->FillH1(3, kinE/keV);   // KinE
    man->FillH1(id, edep/keV);  // по процесах
    man->FillH2(7, kinE/keV, edep/keV);

    // 4) Лог у файл
    fout << name << "\t" << edep/keV << "\n";

    // 5) Передаємо у EventAction
    // знімаємо const-кваліфікатор з базового вказівника
    G4UserEventAction* ua = const_cast<G4UserEventAction*>(
        G4RunManager::GetRunManager()->GetUserEventAction()
        );
    // динамічно приводимо до нашого класу
    if (auto* ea = dynamic_cast<EventAction*>(ua)) {
        ea->AddEdep(edep);
    }

    return true;
  }

private:
  std::ofstream fout;
};

// 3) Геометрія
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    auto* vac  = nist->FindOrBuildMaterial("G4_Galactic");
    auto* NaI  = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    auto* W    = nist->FindOrBuildMaterial("G4_W");

    // World
    G4double R=100*cm;
    auto* solidW = new G4Box("W",R,R,R);
    auto* logicW = new G4LogicalVolume(solidW,vac,"World");
    new G4PVPlacement(nullptr,{},{logicW},"World",nullptr,false,0);

    // Source
    auto* src = new G4Tubs("Src",0,0.5*cm,0.5*mm,0,360*deg);
    auto* lvS = new G4LogicalVolume(src,vac,"Src");
    new G4PVPlacement(nullptr,{0,0,0},lvS,"Src",logicW,false,0);

    // Plate
    G4double zP = 4.0*cm+1.0*mm;
    auto* plt = new G4Box("Plate",5*cm,1.5*cm,1*mm);
    auto* lvP = new G4LogicalVolume(plt,W,"Plate");
    new G4PVPlacement(nullptr,{0,0,zP},lvP,"Plate",logicW,false,0);

    // Scintillator
    G4double zS = zP+1.0*cm+6.3*cm/2;
    auto* sc = new G4Box("Scint",6.3*cm/2,6.3*cm/2,6.3*cm/2);
    auto* lvSc = new G4LogicalVolume(sc,NaI,"Scint");
    new G4PVPlacement(nullptr,{0,0,zS},lvSc,"Scint",logicW,false,0);

    // Підключаємо чутливий детектор
    auto* sd = new SpectrumDetector("SpecSD");
    G4SDManager::GetSDMpointer()->AddNewDetector(sd);
    lvSc->SetSensitiveDetector(sd);

    return new G4PVPlacement(nullptr,{},{logicW},"World",nullptr,false,0);
  }
};

// 4) Генератор Co-60 γ
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction() {
    gun = new G4ParticleGun(1);
    gun->SetParticleDefinition(G4Gamma::GammaDefinition());
  }
  ~PrimaryGeneratorAction() override { delete gun; }
  void GeneratePrimaries(G4Event* evt) override {
    G4double E = (G4UniformRand()<0.5?1173:1332)*keV;
    gun->SetParticleEnergy(E);
    gun->SetParticlePosition({0,0,-5*cm});
    gun->SetParticleMomentumDirection({0,0,1});
    gun->GeneratePrimaryVertex(evt);
    auto* man = G4AnalysisManager::Instance();
    man->FillH1(0, E/keV);  // SourceGamma
  }
private:
  G4ParticleGun* gun;
};

// 5) PhysicsList
class PhysicsList : public G4VModularPhysicsList {
public:
  PhysicsList(){
    RegisterPhysics(new G4EmStandardPhysics());
    RegisterPhysics(new G4DecayPhysics());
    RegisterPhysics(new G4RadioactiveDecayPhysics());
  }
  virtual ~PhysicsList() {}
  virtual void ConstructParticle() override {
      // Стандартні частинки, необхідні для розпаду
      G4Gamma::GammaDefinition();
      G4Electron::ElectronDefinition();
      G4Positron::PositronDefinition();
      G4Proton::ProtonDefinition();
      G4AntiProton::AntiProtonDefinition();
      // Включення іонів
      G4GenericIon::GenericIonDefinition();

      // Обов'язкове створення елементарних частинок для G4RadioactiveDecayPhysics
      G4MuonPlus::MuonPlusDefinition();
      G4MuonMinus::MuonMinusDefinition();
      G4PionPlus::PionPlusDefinition();
      G4PionMinus::PionMinusDefinition();
      G4KaonPlus::KaonPlusDefinition();
      G4KaonMinus::KaonMinusDefinition();
  }
};

// 6) RunAction — усі гістограми
class RunAction : public G4UserRunAction {
public:
  RunAction() {
    std::filesystem::create_directories("results");
    auto* man = G4AnalysisManager::Instance();
    man->SetVerboseLevel(1);
    man->SetDefaultFileType("csv");
    man->SetFileName("results/spectrum");

    // 0–3: SourceGamma, DetectedGamma, EdepFine, KinE
    man->CreateH1("SourceGamma","γ from Co60",200,0.,2000.); // ID=0
    man->CreateH1("DetectedGamma","γ total Edep",200,0.,2000.); //1
    man->CreateH1("EdepFine","Edep fine [keV]",200,0.,1.0); //2
    man->CreateH1("KinE","KinE [keV]",200,0.,2000.); //3

    // 4–6: Edep per process
    man->CreateH1("Photoelectric","Edep photoeffect",200,0.,2000.); //4
    man->CreateH1("Compton","Edep Compton",200,0.,2000.);          //5
    man->CreateH1("PairCreation","Edep pair",200,0.,2000.);        //6

    // 7: 2D — KinE vs Edep
    man->CreateH2("KinE_vs_Edep","KinE vs Edep",200,0.,2000,200,0.,2000); //ID=7

    man->OpenFile();
  }
  void BeginOfRunAction(const G4Run*) override {
    G4AnalysisManager::Instance()->Reset();
  }
  void EndOfRunAction(const G4Run* run) override {
    auto* man = G4AnalysisManager::Instance();
    man->Write();
    man->CloseFile(false);
    G4cout << "Total events: " << run->GetNumberOfEvent() << G4endl;
  }
};

// 7) ActionInitialization
class ActionInitialization : public G4VUserActionInitialization {
public:
  void Build() const override {
    SetUserAction(new PrimaryGeneratorAction());
    SetUserAction(new EventAction());
    SetUserAction(new RunAction());
  }
};

//} // namespace ENX07

// ===== Інтеграція з Qt-мостом або main() =====
extern "C" void InitUserClasses(G4RunManager* run) {
    run->SetUserInitialization(new DetectorConstruction());
    run->SetUserInitialization(new PhysicsList());
    run->SetUserInitialization(new ActionInitialization());
    G4cout << "[ENX07] User classes set\n";
}
namespace Enterprise {
    int StartFromBridge(int argc, char** argv);  // або без namespace, якщо нема
}
int main(int argc, char** argv) {
    return Enterprise::StartFromBridge(argc, argv);
}
