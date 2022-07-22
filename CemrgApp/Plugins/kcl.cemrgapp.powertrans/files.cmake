set(SRC_CPP_FILES

)

set(INTERNAL_CPP_FILES
  kcl_cemrgapp_powertrans_Activator.cpp
  powertransView.cpp
)

set(UI_FILES
  src/internal/powertransViewControls.ui
  src/internal/powertransViewUIRibSpacing.ui
  src/internal/powertransViewUIAhaInput.ui
)

set(MOC_H_FILES
  src/internal/kcl_cemrgapp_powertrans_Activator.h
  src/internal/powertransView.h
)

# list of resource files which can be used by the plug-in
# system without loading the plug-ins shared library,
# for example the icon used in the menu and tabs for the
# plug-in views in the workbench
set(CACHED_RESOURCE_FILES
  resources/icon.xpm
  plugin.xml
)

# list of Qt .qrc files which contain additional resources
# specific to this plugin
set(QRC_FILES

)

set(CPP_FILE