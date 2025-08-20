#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4UserEventAction.hh"
#include "G4VisAttributes.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"

#include "G4VUserActionInitialization.hh"
#include "G4UserRunAction.hh"

#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4GenericIon.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh" // Для G4Proton, G4AntiProton, G4GenericIon

// Включаємо явну ініціалізацію необхідних частинок
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"

#include <windows.h>
#include <fstream>
#include <vector>

#define NUM_CHANNELS 1024

namespace ENX05 {
    	// ------------------- EventAction -------------------
    class EventAction : public G4UserEventAction {
    public:
        void BeginOfEventAction(const G4Event* evt) override {
            edepTotal = 0.0;
            evtID = evt->GetEventID();
        }

        void EndOfEventAction(const G4Event*) override {
            if (edepTotal > 0.)
                G4cout << "Event #" << evtID
                << " — Edep = " << edepTotal / keV << " keV" << G4endl;
        }
        void AddEdep(double edep) { edepTotal += edep; }
    private:
        G4double edepTotal = 0.0;
        G4int evtID = 0;
    };
    // SpectrumDetector: Основний детектор для накопичення спектру та фону.
    class SpectrumDetector : public G4VSensitiveDetector {
    public:
        SpectrumDetector(const G4String& name)
            : G4VSensitiveDetector(name),
            photonCount(0), totalPhotonEnergy(0.), totalEnergyLoss(0.) {
            photonEnergyHist.resize(NUM_CHANNELS, 0);
            energyLossHist.resize(NUM_CHANNELS, 0);
        }

        G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
            G4ParticleDefinition* particleDef = step->GetTrack()->GetDefinition();
            if (particleDef == G4Gamma::GammaDefinition()) {
                double E = step->GetPreStepPoint()->GetKineticEnergy();
                photonCount++;
                totalPhotonEnergy += E;

                int binPhoton = std::clamp(
                    static_cast<int>((E / (2.0 * MeV)) * NUM_CHANNELS),
                    0, NUM_CHANNELS - 1);
                photonEnergyHist[binPhoton]++;

                double r = G4UniformRand();
                double lost = 0.;
                if (r < 0.33) lost = E;
                else if (r < 0.66) lost = E * (0.2 + 0.6 * G4UniformRand());
                else lost = E;

                totalEnergyLoss += lost;

                int binLoss = std::clamp(
                    static_cast<int>((lost / (2.0 * MeV)) * NUM_CHANNELS),
                    0, NUM_CHANNELS - 1);
                energyLossHist[binLoss]++;
            }
            return true;
        }
        void SaveResults() {
            std::ofstream fout1("photon_energy.txt"), fout2("energy_loss.txt");
            for (int i = 0; i < NUM_CHANNELS; ++i) {
                fout1 << photonEnergyHist[i] << "\n";
                fout2 << energyLossHist[i] << "\n";
            }
        }
        void PrintPhotonHistogramAsBar() {
            G4cout << "\nСпектр вхідних енергій (photon_energy.txt):\n";

            // Знайдемо максимальну кількість — для масштабування
            int maxCount = *std::max_element(photonEnergyHist.begin(), photonEnergyHist.end());
            if (maxCount == 0) {
                G4cout << "ℹ️ Жодного фотона не було зафіксовано.\n";
                return;
            }

            for (int i = 0; i < NUM_CHANNELS; ++i) {
                int count = photonEnergyHist[i];
                if (count == 0) continue;

                double energy_keV = (2.0 * i / NUM_CHANNELS) * 1000.0;
                int barLen = static_cast<int>((40.0 * count) / maxCount);  // Масштаб до 40 символів

                G4cout << std::setw(5) << std::fixed << std::setprecision(0)
                    << energy_keV << " keV | "
                    << std::string(barLen, '█') << " "
                    << count << G4endl;
            }
        }
        void PrintEnergyLossHistogramAsBar() {
            G4cout << "=========================================" << G4endl;

            G4cout << "\nСпектр втрат енергії (energy_loss.txt):\n";

            int maxCount = *std::max_element(energyLossHist.begin(), energyLossHist.end());
            if (maxCount == 0) {
                G4cout << "ℹ️ Втрат енергії не зафіксовано.\n";
                return;
            }

            for (int i = 0; i < NUM_CHANNELS; ++i) {
                int count = energyLossHist[i];
                if (count == 0) continue;

                double E = (2.0 * i / NUM_CHANNELS) * 1000.0;  // keV
                int barLen = static_cast<int>((40.0 * count) / maxCount);

                G4cout << std::setw(5) << std::fixed << std::setprecision(0)
                    << E << " keV | "
                    << std::string(barLen, '█') << " "
                    << count << G4endl;
            }

            G4cout << "=========================================" << G4endl;
        }

    private:
        std::vector<int> photonEnergyHist, energyLossHist;
        int photonCount;
        G4double totalPhotonEnergy, totalEnergyLoss;
    };


    // PhysicsList: Фізична модель, що включає електромагнітні процеси та розпад.
    class PhysicsList : public G4VModularPhysicsList {
    public:
        PhysicsList() {
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
    // DetectorConstruction: Геометрія сцинтиляційного детектора.
    class DetectorConstruction : public G4VUserDetectorConstruction {
    public:
        DetectorConstruction() {}
        virtual ~DetectorConstruction() {}

        virtual G4VPhysicalVolume* Construct() override {
            G4NistManager* nist = G4NistManager::Instance();
            G4Material* NaI = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");
            G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");

            G4double worldSize = 100 * cm;
            G4Box* worldS = new G4Box("World", worldSize / 2., worldSize / 2., worldSize / 2.);
            G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, vacuum, "WorldLV");
            G4VPhysicalVolume* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV,
                "WorldPV", nullptr, false, 0, true);

            G4double detSize = 6.3 * cm;
            G4Box* scintS = new G4Box("Scintillator", detSize / 2., detSize / 2., detSize / 2.);
            G4ThreeVector detPosition(0, 0, 3.15 * cm);
            G4LogicalVolume* scintLV = new G4LogicalVolume(scintS, NaI, "ScintillatorLV");
            new G4PVPlacement(nullptr, detPosition, scintLV, "Scintillator", worldLV, false, 0, true);

            SpectrumDetector* mySD = new SpectrumDetector("SpectrumSD");
            G4SDManager::GetSDMpointer()->AddNewDetector(mySD);
            scintLV->SetSensitiveDetector(mySD);
            // Створити візуальні атрибути
            auto* scintAttr = new G4VisAttributes(G4Colour::Red());
                scintAttr->SetForceWireframe(true);
                scintAttr->SetVisibility(true);
                scintLV->SetVisAttributes(scintAttr);  
            auto* worldAttr = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
                worldAttr->SetVisibility(false);
                worldLV->SetVisAttributes(worldAttr); 
            return worldPV;
        }
    };
    // PrimaryGeneratorAction: Джерело частинок Co‑60 (дві гамма-лінії).
    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    public:
        PrimaryGeneratorAction() {
            particleGun = new G4ParticleGun(1);
            particleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
        }
        virtual ~PrimaryGeneratorAction() { delete particleGun; }
        virtual void GeneratePrimaries(G4Event* anEvent) override {
            G4ParticleDefinition* gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
            if (!gamma) {
                G4cerr << "ERROR: gamma particle not found!" << G4endl;
                return;
            }
            particleGun->SetParticleDefinition(gamma);
            GenerateGamma(anEvent, 1173 * keV);
            GenerateGamma(anEvent, 1332 * keV);
        }
        void GenerateGamma(G4Event* anEvent, G4double energy) {
            particleGun->SetParticleEnergy(energy);
            particleGun->SetParticlePosition(G4ThreeVector(0, 0, -5 * cm));
            particleGun->GeneratePrimaryVertex(anEvent);
        }
    private:
        G4ParticleGun* particleGun;
    };
    // RunAction: Запис спектру та фону по завершенні симуляції.
    class RunAction : public G4UserRunAction {
    public:
        virtual void EndOfRunAction(const G4Run*) override {
            SpectrumDetector* specDet =
                dynamic_cast<SpectrumDetector*>(G4SDManager::GetSDMpointer()->FindSensitiveDetector("SpectrumSD"));
            if (auto* det = dynamic_cast<SpectrumDetector*>(
                G4SDManager::GetSDMpointer()->FindSensitiveDetector("SpectrumSD"))) {
                specDet->SaveResults();
                specDet->PrintPhotonHistogramAsBar();      // спектр вхідних фотонів
                specDet->PrintEnergyLossHistogramAsBar();  // спектр втрат

            }
            else
        G4cout << "SpectrumDetector (SpectrumSD) not found!" << G4endl;
        }
    };
    // ActionInitialization: Реєстрація дій (первинний генератор + RunAction).
    class ActionInitialization : public G4VUserActionInitialization {
    public:
        virtual void Build() const override {
            SetUserAction(new PrimaryGeneratorAction());
            SetUserAction(new RunAction());
			SetUserAction(new EventAction()); // 🔧 ДОДАНО!
        }
    };
} // end namespace ENX05
namespace Enterprise
{
    int StartFromBridgeint(argc, char** argv) {
    auto* runManager = new G4RunManager;
    std::cout << "[ENTER] StartFromBridge()\n";
    runManager->SetUserInitialization(new ENX03::DetectorConstruction());
    runManager->SetUserInitialization(new ENX03::PhysicsList());
    runManager->SetUserInitialization(new ENX03::ActionInitialization());
    runManager->Initialize();
    
    return 0;
    }
}