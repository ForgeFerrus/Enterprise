// ENX04.cc - Сцинтиляційний спектрометр для Co-60 (1.3 kV, 100 с)
// ✔ Джерело: Co-60 (γ 1173 keV, 1332 keV)
// ✔ Детектор: NaI(Tl) сцинтилятор
// ✔ Візуалізація як у стандартному прикладі B1
// ✔ Запис спектру на 1024 канали у spectrum.dat
// ✔ Вивід у термінал
// ----------------------------------------------------------------------------
#include "G4Version.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4UserEventAction.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4VisAttributes.hh"

#include <functional>
#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <windows.h> // для SetConsoleOutput
#include <iomanip> // для std::setw

namespace ENX04
{
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
        std::ofstream file("data/spectrum.dat", std::ios::app);
                << " — Edep = " << edepTotal / keV << " keV" << G4endl;
        }
        void AddEdep(double edep) { edepTotal += edep; }
    private:
        G4double edepTotal = 0.0;
        G4int evtID = 0;
    };
    // ------------------- Сцинтиляційний детектор -------------------
    class ScintillationDetector : public G4VSensitiveDetector {
    public:
        static constexpr int BINS = 1024;
        static constexpr double EMIN = 0.0 * keV;
        static constexpr double EMAX = 1500.0 * keV;

        ScintillationDetector(G4String name)
            : G4VSensitiveDetector(name), hist(BINS, 0), binWidth((EMAX - EMIN) / BINS) {
        }
        // Ініціалізація гістограми
        ~ScintillationDetector() override {
            std::ofstream out("spectrum.dat");
            std::cout << "\n  (keV | подій):\n";
            std::cout << "\n Spectrum scintillation (keV | bins):\n";
            for (int i = 0; i < BINS; ++i) {
                double energy = EMIN + i * binWidth;
                out << std::setw(10) << energy / keV << " "
                    << std::setw(10) << hist[i] << "\n";

                // Виводити тільки непорожні бін або кожен N-ий (щоб не засмічувати термінал)
                if (hist[i] > 0 || (i % 64 == 0 && i > 0))
                    std::cout << std::setw(5) << int(energy / keV)
                    << " keV : " << hist[i] << "\n";
            }
        }
        G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
            // 🔹 Отримуємо edep
            G4double edep = step->GetTotalEnergyDeposit();
            if (edep <= 0.) return false;
            // 🔹 Накопичення в гістограму
            int bin = static_cast<int>((edep - EMIN) / binWidth);
            if (bin >= 0 && bin < BINS)
                hist[bin]++;
            // 🔹 Передача в EventAction (для виводу в реальному часі)
            if (auto* eventAction = const_cast<EventAction*>(
                dynamic_cast<const EventAction*>(
                    G4RunManager::GetRunManager()->GetUserEventAction())))
            {
                eventAction->AddEdep(edep);
            }
            return true;
        }
    private:
        std::vector<int> hist;
        double binWidth;
    };
    // ------------------- Геометрія -------------------
    class DetectorConstruction : public G4VUserDetectorConstruction {
    public:
        G4VPhysicalVolume* Construct() override {
            auto* nist = G4NistManager::Instance();
            auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
            auto* NaI = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");

            auto* worldS = new G4Box("World", 1 * m, 1 * m, 1 * m);
            auto* worldLV = new G4LogicalVolume(worldS, vacuum, "World");
            auto* worldPV = new G4PVPlacement(0, {}, worldLV, "World", nullptr, false, 0);
            auto* worldAttr = new G4VisAttributes(G4Colour::Grey());
            worldAttr->SetVisibility(false); // робимо світ прозорим
            worldLV->SetVisAttributes(worldAttr);

            auto* scintS = new G4Box("Scint", 5 * cm, 5 * cm, 5 * cm);
            auto* scintLV = new G4LogicalVolume(scintS, NaI, "Scintillator");
            new G4PVPlacement(0, {}, scintLV, "Scintillator", worldLV, false, 0);
            auto* sdMan = G4SDManager::GetSDMpointer();
            auto* sd = new ScintillationDetector("ScintSD");
            sdMan->AddNewDetector(sd);
            scintLV->SetSensitiveDetector(sd);
            G4cout << "[DEBUG] DetectorConstruction: creating world volume" << G4endl;
            // Візуалізація сцинтилятора
            auto* scintAttr = new G4VisAttributes(G4Colour::Yellow());
            scintAttr->SetForceWireframe(false);
            scintAttr->SetForceSolid(false);
            scintAttr->SetVisibility(true);
            scintLV->SetVisAttributes(scintAttr);

            return worldPV;
        }
    };
    // ------------------- PrimaryGeneratorAction -------------------
    // Генерує гамма‑промені з Co‑60 (випромінювання: 1173 keV або 1332 keV)
    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction() {
        fParticleGun = new G4ParticleGun(1);

        // Створюємо іон Co‑60 (Z=27, A=60, E=0)
        auto* ion = G4IonTable::GetIonTable()->GetIon(27, 60, 0);
        if (!ion) {
            G4cerr << "[ERROR] Co-60 ion not found!" << G4endl;
            exit(1);
        }

        fParticleGun->SetParticleDefinition(ion);
        fParticleGun->SetParticlePosition({0, 0, -10 * cm});
        fParticleGun->SetParticleMomentumDirection({0, 0, 1});
        fParticleGun->SetParticleEnergy(0.0); // фізика розпаду задає енергії
    }

    void GeneratePrimaries(G4Event* event) override {
        fParticleGun->GeneratePrimaryVertex(event); // запускає радіоактивний розпад
    }

    ~PrimaryGeneratorAction() override {
        delete fParticleGun;
    }

private:
    G4ParticleGun* fParticleGun;
};
    // ------------------- PhysicsList -------------------
    class PhysicsList : public G4VModularPhysicsList {
    public:
        PhysicsList() {
            SetVerboseLevel(1);
            RegisterPhysics(new G4EmStandardPhysics());
            RegisterPhysics(new G4DecayPhysics());
            RegisterPhysics(new G4RadioactiveDecayPhysics());
        }
        virtual ~PhysicsList() {}
    };
    // ------------------- RunAction -------------------
    class RunAction : public G4UserRunAction {
    public:
        RunAction() { G4cout << "RunAction: Initialized." << G4endl; }
        virtual ~RunAction() {}
        virtual void BeginOfRunAction(const G4Run* /*run*/) override {
            G4cout << "RunAction: Begin of run." << G4endl;
        }
        virtual void EndOfRunAction(const G4Run* /*run*/) override {
            G4cout << "RunAction: End of run." << G4endl;
        }
    };
    // ------------------- ActionInitialization -------------------
    class ActionInitialization : public G4VUserActionInitialization {
    public:
        void Build() const override {
            auto* eventAction = new EventAction();
            SetUserAction(new PrimaryGeneratorAction());
            SetUserAction(eventAction);
            SetUserAction(new RunAction()); // опціонально
            // Додаємо детектор
            auto* detector = new ScintillationDetector("ScintSD");
            auto* sdManager = G4SDManager::GetSDMpointer();
            sdManager->AddNewDetector(detector);
        }

        void BuildForMaster() const override {
            SetUserAction(new RunAction());
        }
    };
} // namespace ENX04
    // ------------------- Функція запуску симуляції -------------------
    int main() {
    SetConsoleOutputCP(CP_UTF8);
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    auto* runManager = new G4RunManager;
        runManager->SetUserInitialization(new ENX04::DetectorConstruction());
        runManager->SetUserInitialization(new ENX04::PhysicsList());
        runManager->SetUserInitialization(new ENX04::ActionInitialization());
        runManager->Initialize();
    auto* UIm = G4UImanager::GetUIpointer();
    G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();
    // ... запуск візуального інтерфейсу
        ui->SessionStart();
        G4cout << "Simulation completed successfully." << G4endl;
        delete ui;
        delete visManager;
        delete runManager;
        return 0;  
    }

