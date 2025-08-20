// ENX03.cc — інтегрована симуляція Geant4, одиночний файл
// без МТ режиму, з підтримкою візуалізації та користувацьких дій.
// гранична межа 10 000 подій на запуск. до 100 000 для 1М необхідно МТ
// 
// // Цей файл містить всі необхідні включення, визначення класів та функцій
// для створення простої симуляції Geant4, яка включає в себе
// детектор, фізику частинок, генерацію первинних частинок та візуалізацію.
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
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
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
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
#include "G4ionIonisation.hh"
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
// ——————————————————————————————————————————
#include "G4Run.hh"
#include <atomic>
#include <corecrt_math_defines.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <filesystem>
// ——————————————————————————————————————————
namespace ENX03M 
{
    inline double BIN_WIDTH = 0.0;
    inline double userEnergy = 0.0;
    inline int NUM_BINS = 0;
	double maxEnergy = 10.0; // Максимальна енергія для гістограми (у МеВ)
    void SetParameters(double width, double energy, int bins) {
		BIN_WIDTH = width; // Ширина біну для гістограми
		NUM_BINS = bins; // Кількість бінів для гістограми
		userEnergy = energy; // Енергія користувача для гістограми
        maxEnergy = energy; // Максимальна енергія для гістограми
	}
}
// Sensitive Detector
class MySensitiveDetector : public G4VSensitiveDetector {
public:
	static MySensitiveDetector* gEntrySD;
	static MySensitiveDetector* gExitSD;
    std::map<std::string, std::vector<int>> spectra_;

    MySensitiveDetector(const G4String& name) : G4VSensitiveDetector(name) {
        hist_.resize(ENX03M::NUM_BINS, 0);
    }
    G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
        auto pt = step->GetPreStepPoint();
        // 🔻 Знімаємо фільтр — фіксуємо ВСІ частинки
        // (можна залишити коментар для себе)
        // if (pt->GetStepStatus() != fGeomBoundary) return false;

        double E = pt->GetKineticEnergy() / MeV;
        int bin = std::clamp(int(E / ENX03M::maxEnergy * ENX03M::NUM_BINS), 0, ENX03M::NUM_BINS - 1);
        hist_[bin]++;
        const auto& pname = step->GetTrack()->GetParticleDefinition()->GetParticleName();
		// Записуємо енергію та назву частинки у лог-файл
        spectra_[pname].resize(ENX03M::NUM_BINS);
        spectra_[pname][bin]++;
        logFile << pname << " " << E << " MeV\n";
        return true;
    }
    void Save(const std::string& id) {
        std::filesystem::create_directories("results/spectrum");

        // Запис загальної кількості частинок (не по типах)
        {
            std::ofstream out("results/bin_" + id + ".txt");
            for (int count : hist_)
                out << count << "\n";
        }
    }
    void EndOfEvent(G4HCofThisEvent*) override {
        // якщо треба — можна закрити файл після кожного event'у
        logFile.flush();
    }
    const std::vector<int>& Hist() const { return hist_; }
    const std::map<std::string, std::vector<int>>& Spectra() const { return spectra_; }
    private:
	// Історія частинок у бінарному вигляді
	// (для зберігання кількості частинок у кожному біні)
    std::vector<int> hist_;
	std::ofstream logFile{ "hits.txt", std::ios::app }; // Лог-файл для запису частинок
};

MySensitiveDetector* MySensitiveDetector::gEntrySD = nullptr;
MySensitiveDetector* MySensitiveDetector::gExitSD = nullptr;
// ——————————————————————————————————————————
// DetectorConstruction
class MyDetectorConstruction : public G4VUserDetectorConstruction {
public:
    MyDetectorConstruction() : fEntryLV(nullptr), fExitLV(nullptr) {
    }
    G4VPhysicalVolume* Construct() {
        // 1) Матеріали
        auto* nist = G4NistManager::Instance();
        auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
        auto* tungsten = nist->FindOrBuildMaterial("G4_W");

        // 2) Світ
        G4double worldSize = 1.0 * m;
        auto* worldSolid = new G4Box("World", worldSize / 2, worldSize / 2, worldSize / 2);
        auto* worldLV = new G4LogicalVolume(worldSolid, vacuum, "WorldLV");
        auto* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(),
            worldLV, "WorldPV",
            nullptr, false, 0, true);

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
            new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -3 * cm), tgtLV, "TargetPV", worldLV, false, 0);
            tgtLV->SetVisAttributes(wireOnly(G4Colour(0.4, 0.4, 0.4))); // темно-сірий
        }

        /// Cylinder: R=50mm, L=140mm
            auto* cylSolid = new G4Tubs("Cylinder", 0, 50 * mm, 70 * mm, 0, 360 * deg);
            auto* cylLV = new G4LogicalVolume(cylSolid, vacuum, "CylLV");
            new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 7 * cm), cylLV, "CylPV", worldLV, false, 0, true);
            cylLV->SetVisAttributes(wireOnly(G4Colour(0.7, 0.7, 0.7))); // світло-сірий
        
            // Слайси: ±70 мм + додаткове зміщення → уникнення Overlap
            auto* slice = new G4Tubs("Slice", 0, 50 * mm, 0.5 * mm, 0, 360 * deg);

            G4double offset = 0.6 * mm;
            G4double cylZ = 7 * cm;
            G4double entryZ = cylZ - 70 * mm - offset;
            G4double exitZ = cylZ + 70 * mm + offset;

            fEntryLV = new G4LogicalVolume(slice, vacuum, "EntryLV");
            new G4PVPlacement(nullptr, { 0,0,entryZ }, fEntryLV, "EntryPV", worldLV, false, 0, true);
            fEntryLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 0.5, 0.0))); // зелений

            fExitLV = new G4LogicalVolume(slice, vacuum, "ExitLV");
            new G4PVPlacement(nullptr, { 0,0,exitZ }, fExitLV, "ExitPV", worldLV, false, 0, true);
            fExitLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.5, 0.0, 0.0))); // червоний

            return worldPV;
    }
    void ConstructSDandField() override {
        auto* SDM = G4SDManager::GetSDMpointer();
        auto* entrySD = new MySensitiveDetector("EntrySD");
        auto* exitSD = new MySensitiveDetector("ExitSD");

        SDM->AddNewDetector(entrySD);
        SDM->AddNewDetector(exitSD);

        fEntryLV->SetSensitiveDetector(entrySD);
        fExitLV->SetSensitiveDetector(exitSD);

        MySensitiveDetector::gEntrySD = entrySD;
        MySensitiveDetector::gExitSD = exitSD;
    }
    private:
        G4LogicalVolume* fEntryLV;
        G4LogicalVolume* fExitLV;
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
        SetVerboseLevel(1); // рівень виводу інформації
        G4IonConstructor().ConstructParticle(); // конструктор іонів
	}
    void ConstructParticle() override {
        // Створення частинок
        G4Electron::ElectronDefinition();
        G4Gamma::GammaDefinition();
        G4Positron::PositronDefinition();  // ← ось вона — ключ до тиші!
		G4Proton::Proton();
		// Створення нейтронів та альфа-частинок
        G4Neutron::Neutron();
        G4Alpha::Alpha();
		// Створення всіх іонів

    }
    void ConstructProcess() override {
        AddTransportation(); // 🔑 обов’язково перший

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

        // альфа
        auto* alpha = G4Alpha::Alpha();
        helper->RegisterProcess(new G4hMultipleScattering(), alpha);
        helper->RegisterProcess(new G4ionIonisation(), alpha);

        // протон
        auto* p = G4Proton::Proton();
        helper->RegisterProcess(new G4hMultipleScattering(), p);
        helper->RegisterProcess(new G4ionIonisation(), p);

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
        gun_->SetParticleEnergy(ENX03M::userEnergy * MeV);
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
class RunAction : public G4UserRunAction {
public:
    std::atomic<int> counter{ 0 };
    void Increment() {
        counter++;
    }
    void EndOfRunAction(const G4Run*) override {
        using SD = MySensitiveDetector;
        auto* entrySD = SD::gEntrySD;
        auto* exitSD = SD::gExitSD;

        // 🔹 Генерація energy_bins.txt з комами (для Excel)
        std::vector<std::string> binStrs;
        {
            std::ofstream out("energy_bins.txt");
            const double step = ENX03M::maxEnergy / ENX03M::NUM_BINS; // ≈ 0.016992 MeV
            for (int i = 0; i < ENX03M::NUM_BINS; ++i) {
                double E = i * step;
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(2) << E;
                std::string E_str = oss.str();
                std::replace(E_str.begin(), E_str.end(), '.', ',');
                out << E_str << "\n";
                binStrs.push_back(E_str);
            }
        }
        // 🔹 Запис спектра (вихід: кількість на бін)
        auto saveSpectrum = [&](const std::string& fname, const std::vector<int>& vec) {
            std::ofstream out("results/" + fname + ".txt");
            for (int i = 0; i < binStrs.size(); ++i) {
                int count = (i < vec.size()) ? vec[i] : 0;
                out << count << "\n";
            }
        };
        // 🔹 Запис спектрів по частинках
        auto saveSpectra = [&](const std::map<std::string, std::vector<int>>& all) {
            std::filesystem::create_directories("results/spectrum");
            for (const auto& [pname, vec] : all) {
                if (vec.empty()) continue;
                std::ofstream out("results/spectrum/" + pname + ".txt");
                for (int i = 0; i < binStrs.size(); ++i) {
                    int count = (i < vec.size()) ? vec[i] : 0;
                    out << count << "\n";
                }
            }
        };
        if (entrySD) entrySD->Save("entry");
        if (exitSD)  exitSD->Save("exit");
        if (entrySD) saveSpectra(entrySD->Spectra());
        if (entrySD) {
            const auto& spectra = entrySD->Spectra();
            for (const auto& [pname, spec] : spectra) {
                G4cout << "[Run] " << pname << ": " << spec.size()
                    << " bins, total = " << std::accumulate(spec.begin(), spec.end(), 0)
                    << "\n";
            }
        }
        G4cout << "[Run] Spectrum + bins saved.\n";

    }
};
// ── EventAction
class EventAction : public G4UserEventAction {
public:
    std::chrono::steady_clock::time_point lastPrint;
    EventAction(RunAction* runAction) : fRun(runAction) {}

    void BeginOfEventAction(const G4Event* /*evt*/) override {
        // Можна додати логіку ініціалізації події, якщо потрібно
    }

    void EndOfEventAction(const G4Event* /*event*/) override {
        auto now = std::chrono::steady_clock::now();
        if ((now - lastPrint) > std::chrono::seconds(2)) {
            lastPrint = now;

            int seen = fRun->counter.load();
            int total = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed();
            int percent = static_cast<int>((100.0 * seen) / total);

            G4cout << "[Progress] " << percent << "% (" << seen << "/" << total << ")\n";
        }
    }
    private:
    RunAction* fRun;
};
// ── ActionInitialization
class ActionInitialization : public G4VUserActionInitialization {
public:
    void Build() const override 
    {
        SetUserAction(new PrimaryGeneratorAction());

        auto* runAction = new RunAction();
        SetUserAction(runAction);

        auto* eventAction = new EventAction(runAction);
        SetUserAction(eventAction);
    }
};

// ── main()
int main(int argc, char** argv) {
    std::cout << "Enter BIN_WIDTH [MeV]: ";
    std::cin >> ENX03M::BIN_WIDTH;
    ENX03M::NUM_BINS = int(10.0 / ENX03M::BIN_WIDTH);

    std::cout << "Enter primary energy [MeV]: ";
    std::cin >> ENX03M::userEnergy;

    auto* run = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    run->SetUserInitialization(new MyDetectorConstruction());
    run->SetUserInitialization(new MyPhysicsList());
    run->SetUserInitialization(new ActionInitialization());
    run->Initialize();
    // 1) Створюємо й ініціалізуємо візуалізатор
    auto* visManager = new G4VisExecutive;
    visManager->Initialize();

    // 2) Після цього виконуємо макро для візуалізації
    G4UImanager* UIm = G4UImanager::GetUIpointer();
    UIm->ApplyCommand("/control/execute init_vis.mac");

    // 3) Потім стартуємо UI
    auto* ui = new G4UIExecutive(argc, argv);
    ui->SessionStart();
    std::cout << "=== Simulation started ===\n";

    delete ui;
    std::cout << "=== Simulation finished ===\n";
    delete visManager;
    delete run;
    return 0;
}
