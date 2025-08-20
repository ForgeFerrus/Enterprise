
// ENX09МТ.cc Simulation Enterprise project Version 1.0
// Цей файл містить основний код для запуску симуляції Geant4 з підтримкою багатопоточності (MT).
// Він включає налаштування менеджера запуску, фізичних процесів, чутливих детекторів та візуалізації.
// Він також містить налаштування для гістограмування енергії частинок та запису результатів у файл.
// Модель: Електронний пучок 8.7 MeV, вольфрамова пластинка як мішень,
//       чутливий детектор для фіксації спектральних даних.
// Результати записуються у файли:
//    - spectrum.dat  (енергетичний спектр)
//    - angle.dat     (кутовий розподіл)
// ----------------------------------------------------------------------------
#include "G4MTRunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4UserRunAction.hh"
#include "G4MultiRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserStackingAction.hh"
#include "G4UserTrackingAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicsListHelper.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4IonConstructor.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4HadronElasticProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronFissionProcess.hh"
#include "G4ProcessManager.hh"
#include "G4CascadeInterface.hh"
#include "G4LeptonConstructor.hh"
// Моделі та XS-таблиці
#include "G4HadronElastic.hh"
#include "G4NeutronElasticXS.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4LFission.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleChange.hh"
#include "G4ProcessType.hh"
// Параметри користувача
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Threading.hh"
#include <map>
#include <sstream>
#include <chrono>
#include <atomic>
#include <string>
#include <iostream>
#include <iomanip> // для std::setw
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <filesystem>
#include <corecrt_math_defines.h>

// ——————————————————————————————————————————
namespace ENX09
{
    inline double BIN_WIDTH = 0.0;
    inline double userEnergy = 0.0;
    inline int NUM_BINS = 0;
    inline double maxEnergy = 0.0; // Максимальна енергія для гістограми (у МеВ)
    void SetParameters(double width, double energy, int bins) {
        BIN_WIDTH = width; // Ширина біну для гістограми
        NUM_BINS = bins; // Кількість бінів для гістограми
        userEnergy = energy; // Енергія користувача для гістограми
    }
// Sensitive Detector
// Базовий SD з фільтрацією за типом частинки
class MySensitiveDetector : public G4VSensitiveDetector {
public:
    static MySensitiveDetector* gammaEntrySD;
    static MySensitiveDetector* neutronExitSD;
    static MySensitiveDetector* resultSD;
    enum class Type { Gamma, Neutron, Any };
    MySensitiveDetector(const G4String& name, Type type)
        : G4VSensitiveDetector(name), type_(type), hist_(ENX09::NUM_BINS, 0) {}

    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
        auto pt = step->GetPreStepPoint();
        if (pt->GetKineticEnergy() <= 0) return false;
        const auto pname = step->GetTrack()->GetParticleDefinition()->GetParticleName();
        // Фільтрація за типом SD
        if (type_ == Type::Gamma && pname != "gamma") return false;
        if (type_ == Type::Neutron && pname != "neutron") return false;
        double E = pt->GetKineticEnergy() / MeV;
        int bin = std::clamp(int(E / ENX09::maxEnergy * ENX09::NUM_BINS), 0, ENX09::NUM_BINS - 1);
        hist_[bin]++;
        spectra_[pname].resize(ENX09::NUM_BINS);
        spectra_[pname][bin]++;
        return true;
    }
    const std::vector<int>& Hist() const { return hist_; }
    const std::map<std::string, std::vector<int>>& Spectra() const { return spectra_; }
private:
    Type type_;
    std::vector<int> hist_;
    std::map<std::string, std::vector<int>> spectra_;
};

// Глобальний контейнер для всіх SD
inline std::vector<MySensitiveDetector*> allSDs;

// —————————————————————————————————————————————————————————————————
// Оголошення та ініціалізація статичних членів:
MySensitiveDetector* MySensitiveDetector::gammaEntrySD = nullptr;
MySensitiveDetector* MySensitiveDetector::neutronExitSD = nullptr;
MySensitiveDetector* MySensitiveDetector::resultSD = nullptr;
// ——————————————————————————————————————————
// DetectorConstruction
class MyDetectorConstruction : public G4VUserDetectorConstruction {
public:
    MyDetectorConstruction() : fGammaEntryLV(nullptr), fNeutronExitLV(nullptr), fResultLV(nullptr) {}
    G4VPhysicalVolume* Construct() {
        // 1) Матеріали
        auto* nist = G4NistManager::Instance();
        auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
        auto* tungsten = nist->FindOrBuildMaterial("G4_W");
        auto* beryllium = nist->FindOrBuildMaterial("G4_Be");
        auto* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

        // 2) Світ
        G4double worldSize = 1.0 * m;
        auto* worldSolid = new G4Box("World", worldSize / 2, worldSize / 2, worldSize / 2);
        auto* worldLV = new G4LogicalVolume(worldSolid, vacuum, "WorldLV");
        auto* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "WorldPV", nullptr, false, 0, true);

        // 3) Lambda для чистого wireframe
        auto wireOnly = [](G4Colour c) {
            auto* vis = new G4VisAttributes(c);
            vis->SetForceWireframe(true);
            vis->SetForceSolid(false);
            vis->SetForceAuxEdgeVisible(true);
            return vis;
        };

        // 4) Джерело — сфера Gun (радіус 1.5 мм) у z = –10 см
        {
            auto* gunSolid = new G4Sphere("GunSolid", 0, 1.5 * mm, 0, 360 * deg, 0, 180 * deg);
            auto* gunLV = new G4LogicalVolume(gunSolid, vacuum, "GunLV");
            new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -10 * cm), gunLV, "GunPV", worldLV, false, 0, true);
            gunLV->SetVisAttributes(wireOnly(G4Colour(0.4, 0.7, 1.0)));
        }

        // 5) Вольфрамова мішень (93×55×2 мм)
        {
            auto* tgtSolid = new G4Box("Target", 93 / 2 * mm, 55 / 2 * mm, 1 * mm);
            auto* tgtLV = new G4LogicalVolume(tgtSolid, tungsten, "TargetLV");
            new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -5 * cm), tgtLV, "TargetPV", worldLV, false, 0);
            tgtLV->SetVisAttributes(wireOnly(G4Colour(0.4, 0.4, 0.4)));
        }

        // 6) Параметри блоків
        G4double blockThick = 7 * cm;
        G4SDManager* sdm = G4SDManager::GetSDMpointer();

        // 7) Блок No1: Be1 + сталевий кожух 1 мм
        G4ThreeVector posBe1(0, 0, 2 * cm + blockThick / 2);
        auto* be1Solid = new G4Tubs("Be1Core", 0, 5 * cm, blockThick / 2, 0, 360 * deg);
        auto* be1LV = new G4LogicalVolume(be1Solid, beryllium, "Be1LV");
        be1LV->SetVisAttributes(wireOnly(G4Colour(0.2, 0.6, 0.8)));
        new G4PVPlacement(nullptr, posBe1, be1LV, "Be1PV", worldLV, false, 0, true);
        // сталевий кожух
        auto* shell1Solid = new G4Tubs("Shell1", 5 * cm, 5.1 * cm, blockThick / 2, 0, 360 * deg);
        auto* shell1LV = new G4LogicalVolume(shell1Solid, steel, "Shell1LV");
        shell1LV->SetVisAttributes(wireOnly(G4Colour(0.6, 0.6, 0.6)));
        new G4PVPlacement(nullptr, posBe1, shell1LV, "Shell1PV", worldLV, false, 0, true);

        // 8) Блок No2: Be2 + сталевий кожух 2 мм
        G4ThreeVector posBe2(0, 0, 2 * cm + blockThick + 1 * cm + blockThick / 2);
        auto* be2Solid = new G4Tubs("Be2Core", 0, 5 * cm, blockThick / 2, 0, 360 * deg);
        auto* be2LV = new G4LogicalVolume(be2Solid, beryllium, "Be2LV");
        be2LV->SetVisAttributes(wireOnly(G4Colour(0.2, 0.6, 0.8)));
        new G4PVPlacement(nullptr, posBe2, be2LV, "Be2PV", worldLV, false, 0, true);
        auto* shell2Solid = new G4Tubs("Shell2", 5 * cm, 5.2 * cm, blockThick / 2, 0, 360 * deg);
        auto* shell2LV = new G4LogicalVolume(shell2Solid, steel, "Shell2LV");
        shell2LV->SetVisAttributes(wireOnly(G4Colour(0.7, 0.7, 0.7)));
        new G4PVPlacement(nullptr, posBe2, shell2LV, "Shell2PV", worldLV, false, 0, true);

        // 9) Додаємо чутливі області:
        // 9.1. GammaEntrySD — тонкий диск перед Be1
        G4double gammaEntryZ = posBe1.z() - blockThick / 2 - 0.5 * mm;
        auto* gammaEntrySolid = new G4Tubs("GammaEntrySD", 0, 5 * cm, 0.25 * mm, 0, 360 * deg);
        fGammaEntryLV = new G4LogicalVolume(gammaEntrySolid, vacuum, "GammaEntryLV");
        fGammaEntryLV->SetVisAttributes(wireOnly(G4Colour(1.0, 1.0, 0.0)));
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, gammaEntryZ), fGammaEntryLV, "GammaEntryPV", worldLV, false, 0, true);

        // 9.2. NeutronExitSD — тонкий диск після Be2
        G4double neutronExitZ = posBe2.z() + blockThick / 2 + 0.5 * mm;
        auto* neutronExitSolid = new G4Tubs("NeutronExitSD", 0, 5 * cm, 0.25 * mm, 0, 360 * deg);
        fNeutronExitLV = new G4LogicalVolume(neutronExitSolid, vacuum, "NeutronExitLV");
        fNeutronExitLV->SetVisAttributes(wireOnly(G4Colour(1.0, 0.0, 0.0)));
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, neutronExitZ), fNeutronExitLV, "NeutronExitPV", worldLV, false, 0, true);

        // 9.3. ResultSD — кінцевий детектор (диск, наприклад, на +30 см)
        G4double resultZ = neutronExitZ + 30 * cm;
        auto* resultSolid = new G4Tubs("ResultSD", 0, 5 * cm, 0.5 * mm, 0, 360 * deg);
        fResultLV = new G4LogicalVolume(resultSolid, vacuum, "ResultLV");
        fResultLV->SetVisAttributes(wireOnly(G4Colour(0.0, 1.0, 0.0)));
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, resultZ), fResultLV, "ResultPV", worldLV, false, 0, true);
        return worldPV;
    }
    void ConstructSDandField() override {
        auto* SDM = G4SDManager::GetSDMpointer();
        // GammaEntrySD: тільки gamma
        auto* gammaEntrySD = new MySensitiveDetector("GammaEntrySD", MySensitiveDetector::Type::Gamma);
        SDM->AddNewDetector(gammaEntrySD);
        fGammaEntryLV->SetSensitiveDetector(gammaEntrySD);
        allSDs.push_back(gammaEntrySD);
        // NeutronExitSD: тільки neutron
        auto* neutronExitSD = new MySensitiveDetector("NeutronExitSD", MySensitiveDetector::Type::Neutron);
        SDM->AddNewDetector(neutronExitSD);
        fNeutronExitLV->SetSensitiveDetector(neutronExitSD);
        allSDs.push_back(neutronExitSD);
        // ResultSD: всі частинки
        auto* resultSD = new MySensitiveDetector("ResultSD", MySensitiveDetector::Type::Any);
        SDM->AddNewDetector(resultSD);
        fResultLV->SetSensitiveDetector(resultSD);
        allSDs.push_back(resultSD);
    }
private:
    G4LogicalVolume* fGammaEntryLV;
    G4LogicalVolume* fNeutronExitLV;
    G4LogicalVolume* fResultLV;
};
// ——————————————————————————————————————————
// ── MyPhotoNuclearProcess
// Цей клас реалізує процес фотоядерної взаємодії
class MyPhotoNuclearProcess : public G4VDiscreteProcess {
public:

    G4ParticleChange fParticleChange; // для вторинних частинок
    MyPhotoNuclearProcess() : G4VDiscreteProcess("MyPhotoNuclear") {
        SetProcessType(fHadronic);
        SetProcessSubType(131); // fPhotoNuclear
    }
    G4bool IsApplicable(const G4ParticleDefinition* p) { return p == G4Gamma::Gamma(); }
    G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* cond) override {
        *cond = NotForced;
        return 5 * cm; // груба оцінка, або можна зробити енергозалежною
    }

    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step&) override {
        fParticleChange.Initialize(track);
        if (track.GetKineticEnergy() < 7.0 * MeV)
            return &fParticleChange;

        // Створюємо нейтрон як простий вторинний
        auto* neutron = new G4DynamicParticle(G4Neutron::Neutron(), RandomDirection(), 5.0 * MeV);
        fParticleChange.SetNumberOfSecondaries(1);
        fParticleChange.AddSecondary(neutron);

        return &fParticleChange;
    }

private:
    G4ThreeVector RandomDirection() {
        double theta = acos(2 * G4UniformRand() - 1);
        double phi = 2 * M_PI * G4UniformRand();
        return G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
    }
};
// Фізичний список
class MyPhysicsList : public G4VUserPhysicsList {
public:
    MyPhysicsList() {
        // Встановлюємо параметри фізики
        SetVerboseLevel(0); // рівень виводу інформації
    }
    void ConstructParticle() override {
        // Створення частинок
        G4Electron::ElectronDefinition();
        G4Gamma::GammaDefinition();
        G4GenericIon::GenericIonDefinition();
        G4Positron::PositronDefinition();  // ← ось вона — ключ до тиші!
        G4Proton::ProtonDefinition();
        // Створення нейтронів та альфа-частинок
        G4Neutron::Neutron();
        G4Alpha::Alpha();
    }
    void ConstructProcess() override {
        AddTransportation(); // 🔑 обов’язково перший
        // Реєстрація фізичних процесів
        auto* helper = G4PhysicsListHelper::GetPhysicsListHelper();

        // електрони
        auto* eMinus = G4Electron::Electron();
        helper->RegisterProcess(new G4eMultipleScattering(), eMinus);
        helper->RegisterProcess(new G4eIonisation(), eMinus);
        helper->RegisterProcess(new G4eBremsstrahlung(), eMinus);

        // позитрони
        auto* ePlus = G4Positron::Positron();
        helper->RegisterProcess(new G4eMultipleScattering(), ePlus);
        helper->RegisterProcess(new G4eIonisation(), ePlus);
        helper->RegisterProcess(new G4eBremsstrahlung(), ePlus);

        // гамма
        auto* gamma = G4Gamma::Gamma();
        helper->RegisterProcess(new G4PhotoElectricEffect(), gamma);
        helper->RegisterProcess(new G4ComptonScattering(), gamma);
        helper->RegisterProcess(new G4GammaConversion(), gamma);
        // наш кастомний γ → n процес
        auto* myPhoto = new MyPhotoNuclearProcess(); // з класу вище
        helper->RegisterProcess(myPhoto, gamma);

        // альфа-частинки
        auto* alpha = G4Alpha::Alpha();
        helper->RegisterProcess(new G4hMultipleScattering(), alpha);
        helper->RegisterProcess(new G4hIonisation(), alpha);

        // протон
        auto* p = G4Proton::Proton();
        helper->RegisterProcess(new G4hMultipleScattering(), p);
        helper->RegisterProcess(new G4hIonisation(), p);

        // нейтрон
        auto* n = G4Neutron::Neutron();
        auto* elastic = new G4HadronElasticProcess();
        elastic->RegisterMe(new G4HadronElastic());
        elastic->AddDataSet(new G4NeutronElasticXS());
        helper->RegisterProcess(elastic, n);

        auto* capture = new G4NeutronCaptureProcess();
        capture->RegisterMe(new G4NeutronRadCapture());
        capture->AddDataSet(new G4NeutronCaptureXS());
        helper->RegisterProcess(capture, n);

        auto* fission = new G4NeutronFissionProcess();
        fission->RegisterMe(new G4LFission());
        helper->RegisterProcess(fission, n);
    }

    void SetCuts() override {
        for (const auto& name : { "gamma", "e-", "e+", "neutron", "proton", "alpha" })
            SetCutValue(0.001 * mm, name);
    }
};
// ── Генератор
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction() {
        gun_ = new G4ParticleGun(1);
        gun_->SetParticleDefinition(G4Electron::ElectronDefinition());
        gun_->SetParticleEnergy(ENX09::userEnergy * MeV);
        gun_->SetParticleMomentumDirection({ 0, 0, 1 });
        gun_->SetParticlePosition({ 0, 0, -10 * cm });
    }
    ~PrimaryGeneratorAction() override { delete gun_; }
    void GeneratePrimaries(G4Event* evt) override {
        gun_->GeneratePrimaryVertex(evt);
    }

private:
    G4ParticleGun* gun_;
};
// ── RunAction
class RunAction : public G4UserRunAction {
public:
    std::atomic<int> counter = 0;
    void Increment() { counter++; }
    void EndOfRunAction(const G4Run*) override {
        // Генерація energy_bins.txt з комами (для Excel)
        std::filesystem::create_directories("results");
        std::vector<std::string> binStrs;
        {
            std::ofstream out("results/energy_bins.txt");
            const double step = ENX09::maxEnergy / ENX09::NUM_BINS;
            for (int i = 0; i < ENX09::NUM_BINS; ++i) {
                double E = i * step;
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(2) << E;
                std::string E_str = oss.str();
                std::replace(E_str.begin(), E_str.end(), '.', ',');
                out << E_str << "\n";
                binStrs.push_back(E_str);
            }
        }

        // Автоматичний запис спектрів для всіх SD через глобальний контейнер
        for (auto* mysd : allSDs) {
            if (!mysd) continue;
            std::string sdName = mysd->GetName();
            // Гістограма по всім частинкам
            {
                std::ofstream out("results/" + sdName + ".txt");
                const auto& hist = mysd->Hist();
                for (int i = 0; i < binStrs.size(); ++i) {
                    int count = (i < hist.size()) ? hist[i] : 0;
                    out << count << "\n";
                }
            }
        }
    }
};
// ProgressEventAction — просто прогрес-бар у консолі
class ProgressEventAction : public G4UserEventAction {
public:
    ProgressEventAction() : lastPercent_(-1) {}
    virtual ~ProgressEventAction() {}

    virtual void EndOfEventAction(const G4Event* /*evt*/) override {
        // атомарний лічильник усіх подій зрушуємо вгору
        static std::atomic<int> eventCounter{ 0 };
        int count = ++eventCounter;

        // загальна кількість запланованих подій
        int total = G4RunManager::GetRunManager()
            ->GetNumberOfEventsToBeProcessed();

        // обчислюємо відсоток, захищаючи від ділення на нуль
        int percent = (total > 0) ? (count * 100 / total) : 0;

        // друкуємо тільки коли змінюється відсоток
        if (percent != lastPercent_ && count <= total) {
            lastPercent_ = percent;
            G4cout
                << "[Progress] "
                << percent << "% complete ("
                << count << "/" << total << ")"
                << G4endl;
        }
    }

private:
    int lastPercent_;
};
// ── ActionInitialization
class ActionInitialization : public G4VUserActionInitialization {
public:
    void BuildForMaster() const override {
        SetUserAction(new RunAction()); // тільки RunAction для Master
    }

    void Build() const override {
        // спочатку — твій RunAction, PrimaryGeneratorAction тощо
        SetUserAction(new RunAction());
        SetUserAction(new ProgressEventAction());
        SetUserAction(new PrimaryGeneratorAction());
    }
};
} // namespace ENX09M
// —————————————————————————————————————————— 
// Оголошення Qt-режиму без додаткових файлів
namespace Enterprise {
    int StartFromBridgeQt(int argc, char** argv);
}
// ── Головна функція main() — точка входу
    int main(int argc, char** argv) {
        G4cout << "===================================================" << G4endl;
        G4cout << " ENX09 Simulation Program Version 1.0" << G4endl;
        G4cout << " == Enterprise project Version 9.0 ==" << G4endl;
        G4cout << " Modeling particle interactions using a tungsten target" << G4endl;
        G4cout << " and a sensitive detector to record spectral data." << G4endl;
        G4cout << "===================================================" << G4endl;
        // Вибір режиму (MT чи серійний)
        std::cout << "Select mode: [1] Multi-threaded (fast, no interactive visualization)\n"
                     "             [2] Serial (single-threaded, with interactive visualization)\n"
                     "Enter 1 or 2: ";
        int mode = 1;
        std::cin >> mode;
        bool useMT = (mode == 1);

        // Ввід параметрів
        std::cout << "Enter BIN_WIDTH [MeV]: ";
        std::cin >> ENX09::BIN_WIDTH;

        std::cout << "Enter primary energy [MeV]: ";
        std::cin >> ENX09::userEnergy;

        std::cout << "Enter max energy for histogram [MeV]: ";
        std::cin >> ENX09::maxEnergy;

        // Тільки тепер рахуємо кількість бінів
        ENX09::NUM_BINS = std::max(1, int(std::ceil(ENX09::maxEnergy / ENX09::BIN_WIDTH)));

        std::cout << "Enter number of runs: ";
        long long nEvents; std::cin >> nEvents;

        G4RunManager* runManager = nullptr;
        if (useMT) {
            runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
            runManager->SetNumberOfThreads(10);
        } else {
            runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
        }
        runManager->SetUserInitialization(new ENX09::MyDetectorConstruction());
        runManager->SetUserInitialization(new ENX09::MyPhysicsList());
        runManager->SetUserInitialization(new ENX09::ActionInitialization());
        runManager->Initialize();
        runManager->BeamOn(nEvents);
        G4cout << "==== Simulation completed. ====\n" << G4endl;

        // 1) MT: запускаємо Qt-Bridge (за наявності), не видаляючи runManager
#ifdef USE_QT_BRIDGE
        if (useMT) {
            std::cout << "[INFO] Launching Qt-Bridge visualization...\n";
            int qtResult = 0;
            try {
                qtResult = Enterprise::StartFromBridgeQt(argc, argv);
            }
            catch (...) {
                std::cout << "[WARN] Qt-Bridge visualization failed.\n";
                qtResult = -999;
            }
            std::cout << "[DEBUG] Qt-Bridge return code: " << qtResult << std::endl;
        }
#endif

        // 2) Serial: стандартна Geant4 візуалізація (як було)
        if (!useMT) {
            G4VisManager* visManager = new G4VisExecutive;
            visManager->Initialize();

            G4UIExecutive* ui = new G4UIExecutive(argc, argv);
            G4UImanager* UImanager = G4UImanager::GetUIpointer();
            UImanager->ApplyCommand("/control/execute init_vis.mac");
            ui->SessionStart();
            delete ui;
            delete visManager;
        }

        // Тепер можна звільняти
        delete runManager;

        std::cout << "[INFO] Program finished. Press Enter to exit..." << std::endl;
        std::cin.get();
        return 0;
        }

