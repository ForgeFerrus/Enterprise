#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4RunManagerFactory.hh"
#include "G4VisExecutive.hh"
#include <iostream>
#include <string>

#include <DetectorConstruction.hh>
#include <ActionInitialization.hh>
#include <PhysicsList.hh>

namespace Enterprise
{ 
    using namespace std;
    int main(int argc, char** argv) {
        auto* runManager = new G4RunManager();
        auto* detector = new DetectorConstruction();

        runManager->SetUserInitialization(detector);
        runManager->SetUserInitialization(new PhysicsList());

        // Зчитування енергії частинок з консолі
        G4double particleEnergy;
        detector->ReadParticleEnergy(particleEnergy);

        runManager->SetUserInitialization(new ActionInitialization(particleEnergy));

        auto visManager = new G4VisExecutive();
        visManager->Initialize();

        auto ui = new G4UIExecutive(argc, argv);

        G4UImanager* UImanager = G4UImanager::GetUIpointer();
        UImanager->ApplyCommand("/run/initialize");

        if (ui->IsGUI()) {
            UImanager->ApplyCommand("/control/execute init_vis.mac");
        }
        else {
            UImanager->ApplyCommand("/control/execute init.mac");
        }

        // Введення кількості запусків з консолі
        G4int numberOfRuns;
        cout << "Enter number of runs: ";
        cin >> numberOfRuns;

        // Запуск подій
        std::string runCommand = "/run/beamOn " + std::to_string(numberOfRuns);
        UImanager->ApplyCommand(runCommand);

        ui->SessionStart();

        delete visManager;
        delete runManager;
        return 0;
    }

} // закриває простір імені перевіряємо функцію main для виводу в консоль
int main(int argc, char** argv) {
    return Enterprise::main(argc, argv);
}
