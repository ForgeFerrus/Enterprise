// Сцинтиляційний спектрометр для Co-60 (1.3 kV, 100 с)
// 
// Характеристики (лінійні прототипи для лабораторного спектрометра SBS-40):
// • Джерело: Co-60 (γ 1173 keV та 1332 keV)
// • Детектор: кристал NaI(Tl) розміром 63×63 mm (тобто, 63 mm x 63 mm x 63 mm, задаємо як куб)
// • Чутливий детектор: розбиття даних на 1024 канали з розділенням значень енергії 
//   та кількістю зареєстрованих частинок у окремих файлах (spectrum_energy.dat та spectrum_counts.dat)
// • Візуалізація – як у стандартному прикладі B1
//
// ----------------------------------------------------------------------------
// ENX06.cc - Повний код симуляції для Co‑60m із візуалізацією моделі
// (Сцинтиляційний спектрометр для Co‑60, із використанням збудженого стану Co‑60m)
#pragma once
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"
#include "G4IonTable.hh"
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4EmExtraPhysics.hh"
#include "G4GenericIon.hh"

#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <string>

namespace ENX06 
{
    // ------------------- ScintillationDetector -------------------
    // Реалізовує чутливий детектор, що розбиває накопиченні дані на 1024 канали.
    // При завершенні симуляції дані записуються у файли "energy_grad.dat" та "spectrum_Co60.dat".
    class ScintillationDetector : public G4VSensitiveDetector {
    public:
        static const int NOBINS = 1024;
        static const double HIST_MIN;
        static const double HIST_MAX;

        ScintillationDetector(G4String name) : G4VSensitiveDetector(name) {
            for (int i = 0; i < NOBINS; i++) { histogram[i] = 0; }
        }

        virtual ~ScintillationDetector() {
            std::ofstream energyFile("energy_grad.dat");
            std::ofstream countFile("spectrum_Co60.dat");
            double bin_width = (HIST_MAX - HIST_MIN) / NOBINS;
            for (int i = 0; i < NOBINS; ++i) {
                double energy = i * bin_width + HIST_MIN;
                energyFile << std::setw(10) << energy / keV << std::endl;
                countFile << std::setw(10) << histogram[i] << std::endl;
            }
            energyFile.close();
            countFile.close();
            G4cout << "ScintillationDetector: Spectrum data saved." << G4endl;
        }

        virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
            double edep = step->GetTotalEnergyDeposit();
            if (edep <= 0.) return false;
            double bin_width = (HIST_MAX - HIST_MIN) / NOBINS;
            int index = int((edep - HIST_MIN) / bin_width);
            if (index >= 0 && index < NOBINS)
                histogram[index]++;
            G4cout << "Edep: " << edep / keV << " keV" << G4endl;
            return true;
        }

    private:
        int histogram[NOBINS];
    };

    const double ScintillationDetector::HIST_MIN = 0.0;
    const double ScintillationDetector::HIST_MAX = 1500.0 * keV;

    // ------------------- DetectorConstruction -------------------
    // Визначає геометрію симуляції: створює світовий об'єм і розташовує сцинтилятор.
    // Сцинтилятор представлен у формі куба розміром 63 mm × 63 mm × 63 mm.
    class DetectorConstruction : public G4VUserDetectorConstruction {
    public:
        virtual G4VPhysicalVolume* Construct() override {
            G4NistManager* nist = G4NistManager::Instance();
            G4Material* NaI = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE"); // матеріал сцинтилятора
            G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");

            // Світовий об'єм – куб розміром 100 см
            G4double worldSize = 100 * cm;
            G4Box* worldS = new G4Box("World", worldSize, worldSize, worldSize);
            G4LogicalVolume* worldLV = new G4LogicalVolume(worldS, vacuum, "World");
            G4VPhysicalVolume* worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "World", 0, false, 0);

            // Сцинтилятор – куб розміром 63 mm (половинний розмір 3.15 cm)
            G4Box* scintS = new G4Box("Scintillator", 3.15 * cm, 3.15 * cm, 3.15 * cm);
            G4LogicalVolume* scintLV = new G4LogicalVolume(scintS, NaI, "Scintillator");
            new G4PVPlacement(0, G4ThreeVector(0, 0, 0), scintLV, "Scintillator", worldLV, false, 0);

            // Прив'язка чутливого детектора до сцинтилятора
            G4SDManager* sdManager = G4SDManager::GetSDMpointer();
            ScintillationDetector* mySD = new ScintillationDetector("ScintillationDetector");
            sdManager->AddNewDetector(mySD);
            scintLV->SetSensitiveDetector(mySD);

            return worldPV;
        }
    };

    // ------------------- PrimaryGeneratorAction -------------------
    // Генерує первинні частинки.
    // Отримання визначення іона Co‑60m відкладене до першої генерації (GeneratePrimaries),
    /// що гарантує завершення ініціалізації PhysicsList.
    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    public:
        PrimaryGeneratorAction() : fParticleGun(new G4ParticleGun(1)), fParticleDef(nullptr) {
            G4ParticleDefinition* gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
            fParticleGun->SetParticleDefinition(gamma);
            // Розташування джерела позаду детектора (сцинтилятора)
            fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -5 * cm));
            fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));

        }
        virtual ~PrimaryGeneratorAction() { delete fParticleGun; }
        virtual void GeneratePrimaries(G4Event* anEvent) override {
            fParticleGun->SetParticleEnergy(1173 * keV);
            fParticleGun->SetParticleEnergy(1332 * keV);
            fParticleGun->GeneratePrimaryVertex(anEvent);
            if (!fParticleDef) {
                G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
                G4IonTable* ionTable = particleTable->GetIonTable();
                fParticleDef = ionTable->GetIon(27, 60, 1); // Отримання Co-60m
                if (!fParticleDef) {
                    G4Exception("PrimaryGeneratorAction", "IonError", FatalException, "Co-60m not found!");
                }
                fParticleGun->SetParticleDefinition(fParticleDef);
            }
        }
    private:
        G4ParticleGun* fParticleGun;
        G4ParticleDefinition* fParticleDef;
    };
    // ------------------- PhysicsList -------------------
    // Налаштовує фізичні процеси та забезпечує ініціалізацію GenericIon.
    class PhysicsList : public G4VModularPhysicsList {
    public:
        PhysicsList() {
            SetVerboseLevel(1);
            RegisterPhysics(new G4EmStandardPhysics());         // Електромагнітна фізика для гамма-квантів
            RegisterPhysics(new G4DecayPhysics());              // Фізика розпаду
            RegisterPhysics(new G4RadioactiveDecayPhysics());   // Радіоактивний розпад
            RegisterPhysics(new G4EmExtraPhysics());
            RegisterPhysics(new G4HadronPhysicsQGSP_BIC_AllHP());

        }
        virtual ~PhysicsList() {}
        virtual void ConstructParticle() override {
            G4GenericIon::GenericIonDefinition();
            // Викликаємо стандартну ініціалізацію частинок із модульної фізики,
            // Ініціалізація іона Co-60m
            G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(27, 60, 0);
            if (!ion) {
                G4cout << "PhysicsList: Co-60m ion not found!" << G4endl;
            }
            else {
                G4cout << "PhysicsList: Co-60m ion initialized." << G4endl;
            }
        }
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
    // ------------------- PrimaryGeneratorAction -------------------  
    // Генерує первинні частинки.  
    // Отримання визначення іона Co‑60m відкладене до першої генерації (GeneratePrimaries),  
    // що гарантує завершення ініціалізації PhysicsList.  
    class MyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {  
    public:  
       MyPrimaryGeneratorAction() : fParticleGun(nullptr), fParticleDef(nullptr) {  
           fParticleGun = new G4ParticleGun(1);  
           G4ParticleDefinition* gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");  
           fParticleGun->SetParticleDefinition(gamma);  
           // Розташування джерела позаду детектора (сцинтилятора)  
           fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -5 * cm));  
           fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));  
       }  
       virtual ~MyPrimaryGeneratorAction() { delete fParticleGun; }  
       virtual void GeneratePrimaries(G4Event* anEvent) override {  
           fParticleGun->SetParticleEnergy(1173 * keV);  
           fParticleGun->SetParticleEnergy(1332 * keV);  
           fParticleGun->GeneratePrimaryVertex(anEvent);  
           if (!fParticleDef) {  
               G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();  
               G4IonTable* ionTable = particleTable->GetIonTable();  
               fParticleDef = ionTable->GetIon(27, 60, 1); // Отримання Co-60m  
               if (!fParticleDef) {  
                   G4Exception("PrimaryGeneratorAction", "IonError", FatalException, "Co-60m not found!");  
               }  
               fParticleGun->SetParticleDefinition(fParticleDef);  
           }  
       }  
    private:  
       G4ParticleGun* fParticleGun;  
       G4ParticleDefinition* fParticleDef;  
    };
    // ------------------- RunSimulation -------------------
    // Функція задає параметри симуляції, ініціалізує RunManager, запускає візуалізацію та симуляцію.
    int RunSimulation(int argc, char** argv) {
        double voltage;
        std::cout << "Enter detector voltage (kV): ";
        std::cin >> voltage;

        double exposureTime;
        std::cout << "Enter exposure time (s): ";
        std::cin >> exposureTime;

        double eventTime = 0.01;
        int numEvents = std::round(exposureTime / eventTime);
        std::cout << "\nStarting spectrometer simulation...\n";
        std::cout << "Calculated number of events for beamOn: " << numEvents << std::endl;

        // Ініціалізація RunManager
        G4RunManager* runManager = new G4RunManager();
        runManager->SetUserInitialization(new PhysicsList());
        runManager->SetUserInitialization(new DetectorConstruction());
        runManager->Initialize();

        // Ініціалізація і конфігурація візуалізації
        G4VisManager* visManager = new G4VisExecutive();
        visManager->Initialize();
        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        UImanager->ApplyCommand("/vis/open OGL");
        UImanager->ApplyCommand("/vis/scene/create");
        UImanager->ApplyCommand("/vis/drawVolume");
        UImanager->ApplyCommand("/vis/viewer/set/style surface");
        UImanager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 120 150");

        // Опціональні команди, якщо визначені (перевірте їх, якщо не використовуєте)
        std::string commandVoltage = "/detector/set/voltage " + std::to_string(voltage) + "kV";
        UImanager->ApplyCommand(commandVoltage.c_str());
        std::string commandTime = "/detector/set/exposureTime " + std::to_string(exposureTime) + "s";
        UImanager->ApplyCommand(commandTime.c_str());

        // Запуск симуляції
        UImanager->ApplyCommand("/run/beamOn " + std::to_string(numEvents));
        std::cout << "Data saved. Simulation completed.\n";

        // Запуск інтерактивної UI-сесії
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        ui->SessionStart();

        std::cout << "Press ENTER to exit..." << std::endl;
        std::cin.ignore();
        std::cin.get();

        delete ui;
        delete visManager;
        delete runManager;

        return 0;
    }
} // namespace ENX06
// Точка входу в програму
int main(int argc, char** argv) {
    return ENX06::RunSimulation(argc, argv);
}

