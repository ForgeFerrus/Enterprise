#include <G4RunManagerFactory.hh>
#include <G4MTRunManager.hh>
#include <G4NistManager.hh>
#include <G4Threading.hh>
#include <G4RunManager.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4SystemOfUnits.hh>
#include <G4VisAttributes.hh>
#include <G4VisManager.hh>
#include <G4UImanager.hh>
#include <G4AnalysisManager.hh>
#include <G4VUserActionInitialization.hh>
#include <G4VUserPrimaryGeneratorAction.hh>

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4SystemOfUnits.hh>
#include <G4ParticleGun.hh>
#include <G4Event.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VSensitiveDetector.hh>
#include <G4SDManager.hh>
#include <G4Step.hh>
#include <G4UserRunAction.hh>
#include <G4Run.hh>
#include <G4VModularPhysicsList.hh>
#include <G4DecayPhysics.hh>
#include <G4EmStandardPhysics.hh>
#include <G4HadronPhysicsFTFP_BERT.hh>
#include <G4VisExecutive.hh>
#include <G4UImanager.hh>
#include <fstream>
#include <map>
#include <filesystem>

namespace ENX01 {

    class DetectorSD : public G4VSensitiveDetector {
    public:
    DetectorSD(const G4String& name)
        : G4VSensitiveDetector(name),
        HIST_MAX(100 * MeV),
        HIST_MIN(0 * MeV)
    {
        std::fill(std::begin(histogram), std::end(histogram), 0);
        std::fill(std::begin(histogram_angle), std::end(histogram_angle), 0);
    }

    ~DetectorSD() override {
        G4cout << "[DetectorSD] Writing histogram files\n";
        std::ofstream spec("spectrum.dat");
        double binE = (HIST_MAX - HIST_MIN) / NOBINS;
        for (int i = 0; i < NOBINS; ++i) {
            double energy = i * binE + HIST_MIN;
            spec << std::setw(15) << energy / MeV << " " << std::setw(15) << histogram[i] << "\n";
        }
        spec.close();

        std::ofstream ang("angle.dat");
        double binA = CLHEP::pi / NOBINS;
        for (int i = 0; i < NOBINS; ++i) {
            double theta = i * binA;
            ang << std::setw(15) << theta << " " << std::setw(15) << histogram_angle[i] << "\n";
        }
        ang.close();
    }

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
        auto* track = step->GetTrack();
        if (track->GetDefinition()->GetParticleName() != "e-") return false;

        double E = step->GetPreStepPoint()->GetKineticEnergy();
        int eBin = int((E - HIST_MIN) / ((HIST_MAX - HIST_MIN) / NOBINS));
        if (eBin >= 0 && eBin < NOBINS) histogram[eBin]++;

        G4ThreeVector dir = step->GetPreStepPoint()->GetMomentumDirection();
        double angle = dir.angle(G4ThreeVector(0, 0, 1));
        int aBin = int(angle / (CLHEP::pi / NOBINS));
        if (aBin >= 0 && aBin < NOBINS) histogram_angle[aBin]++;

        G4cout << "[Hit] E = " << E / MeV << " MeV, θ = " << angle << " rad\n";
        track->SetTrackStatus(fStopAndKill);
        return true;
    }

private:
    static constexpr int NOBINS = 1000;
    const double HIST_MAX, HIST_MIN;
    int histogram[NOBINS];
    int histogram_angle[NOBINS];
};

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
        DetectorConstruction() = default;
        ~DetectorConstruction() override = default;
    G4LogicalVolume* GetBremsVolume() const { return fBremsVolume; }
    G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
    auto* solidWorld = new G4Box("World", 15*cm, 15*cm, 10*cm);
    auto* logicWorld = new G4LogicalVolume(solidWorld, vacuum, "World");
    auto* physWorld = new G4PVPlacement(nullptr, {}, logicWorld, "World", nullptr, false, 0);

    // Візуалізація: тільки контур (wireframe)
    auto* vis = new G4VisAttributes(G4Colour::Grey());
    vis->SetForceWireframe(true);
    logicWorld->SetVisAttributes(vis);

    // Target
    auto* W = nist->FindOrBuildMaterial("G4_W");
    auto* solidTarget = new G4Tubs("Target", 0, 0.5*cm, 0.5*mm, 0, 360*deg);
    auto* logicTarget = new G4LogicalVolume(solidTarget, W, "Target");
    new G4PVPlacement(nullptr, {}, logicTarget, "Target", logicWorld, false, 0);
    logicTarget->SetVisAttributes(vis);
    fBremsVolume = logicTarget;

    // Detector
    auto* solidDet = new G4Box("Detector", 5*cm, 5*cm, 2.5*cm);
    auto* logicDet = new G4LogicalVolume(solidDet, W, "Detector");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 2.5*cm), logicDet, "Detector", logicWorld, false, 0);
    logicDet->SetVisAttributes(vis);
    logicDetector = logicDet;

    return physWorld;
  }

  void ConstructSDandField() override {
    auto* sd = new DetectorSD("SensitiveDetector");
    G4SDManager::GetSDMpointer()->AddNewDetector(sd);
    logicDetector->SetSensitiveDetector(sd);
  }

private:
  G4LogicalVolume* fBremsVolume = nullptr;
  G4LogicalVolume* logicDetector = nullptr;
};

class PhysicsList : public G4VModularPhysicsList {
public:
    PhysicsList() {
        RegisterPhysics(new G4EmStandardPhysics()); // Стандартна EM фізика
        // При потребі: додай адронну чи кастомну фізику тут
        // RegisterPhysics(new G4HadronPhysicsQGSP_BERT());
    }

    ~PhysicsList() override = default;

    void ConstructParticle() override {
        G4VModularPhysicsList::ConstructParticle();
    }

    void ConstructProcess() override {
        G4VModularPhysicsList::ConstructProcess();
        // AddTransportation(); // Якщо хочеш ручне керування транспортом
        // Можна додати власні процеси тут
    }
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    G4ParticleGun* gun = nullptr;
public:
    PrimaryGeneratorAction() {
        gun = new G4ParticleGun(1);
        gun->SetParticleDefinition(
            G4ParticleTable::GetParticleTable()->FindParticle("e-"));
        gun->SetParticleEnergy(8.7 * MeV);
        gun->SetParticlePosition({ 0,0,-300 * mm });
        gun->SetParticleMomentumDirection({ 0,0,1 });
    }
    ~PrimaryGeneratorAction() override { delete gun; }
    void GeneratePrimaries(G4Event* evt) override {
        G4cout << "Generating vertex…" << G4endl;
        gun->GeneratePrimaryVertex(evt);
    }
};
class SteppingAction : public G4UserSteppingAction {
public:
    SteppingAction() = default;
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override {
        if (!fBremsVolume) {
            auto* detConstruction = static_cast<const DetectorConstruction*>(
                G4RunManager::GetRunManager()->GetUserDetectorConstruction());
            fBremsVolume = detConstruction->GetBremsVolume();
        }

        auto* currentVolume = step->GetPreStepPoint()->GetTouchableHandle()
            ->GetVolume()->GetLogicalVolume();
        if (currentVolume != fBremsVolume) return;

        auto* secondaries = step->GetSecondaryInCurrentStep();
        if (secondaries->empty()) return;

        double electronEnergy = step->GetPreStepPoint()->GetKineticEnergy();
        auto* analysisManager = G4AnalysisManager::Instance();

        for (const auto* track : *secondaries) {
            if (track->GetDefinition()->GetParticleName() != "gamma") continue;

            double gammaE = track->GetKineticEnergy();
            analysisManager->FillNtupleDColumn(0, 0, gammaE);
            analysisManager->AddNtupleRow(0);

            double relE = gammaE / electronEnergy;
            analysisManager->FillNtupleDColumn(1, 0, relE);
            analysisManager->AddNtupleRow(1);
        }
    }

private:
    G4LogicalVolume* fBremsVolume = nullptr;
};
class RunAction : public G4UserRunAction {
public:
    RunAction() {
        auto* analysisManager = G4AnalysisManager::Instance();
        analysisManager->SetVerboseLevel(1);
        analysisManager->CreateNtuple("AbsGamma", "Gamma from target");
        analysisManager->CreateNtupleDColumn(0, "Energy");
        analysisManager->CreateNtuple("RelGamma", "Gamma/electron ratio");
        analysisManager->CreateNtupleDColumn(1, "RelEnergy");
        analysisManager->FinishNtuple();
    }

    ~RunAction() override {
        delete G4AnalysisManager::Instance();
    }

    void BeginOfRunAction(const G4Run* /*run*/) override {
        auto* analysisManager = G4AnalysisManager::Instance();
        analysisManager->OpenFile("output.root");
        G4cout << "[RunAction] Simulation run started" << G4endl;
    }

    void EndOfRunAction(const G4Run* /*run*/) override {
        auto* analysisManager = G4AnalysisManager::Instance();
        analysisManager->Write();
        analysisManager->CloseFile();
        G4cout << "[RunAction] Run ended, data written to output.root" << G4endl;
    }
};
class ActionInitialization : public G4VUserActionInitialization {
public:
    void Build() const override {
        SetUserAction(new PrimaryGeneratorAction);
        SetUserAction(new RunAction);
        SetUserAction(new SteppingAction);
        // інші дії, якщо потрібно
    }

    void BuildForMaster() const override {
        SetUserAction(new RunAction); // тільки те, що потрібно на master
    }
};

} // namespace ENX01
namespace Enterprise {
    int StartFromBridge(int argc, char** argv);  // або без namespace, якщо нема
}
int main(int argc, char** argv) {
    // Вивід тексту
    std::wcout << L"Hello, world!" << std::endl;
	// Ініціалізація менеджера запуску
    auto* run = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);
    run->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
    run->SetUserInitialization(new ENX01::DetectorConstruction());
    run->SetUserInitialization(new ENX01::PhysicsList());
	run->SetUserInitialization(new ENX01::ActionInitialization());
    run->Initialize();

    delete run;
    return Enterprise::StartFromBridge(argc, argv);
}
