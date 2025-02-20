// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIUsersdIcz263673dIalicedIKit_AnalysisdIFlowdIFlowAnalysis_Helper_cxx_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/Users/cz263673/alice/Kit_Analysis/Flow/./FlowAnalysis_Helper.cxx"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *FlowAnalysis_Helper_Dictionary();
   static void FlowAnalysis_Helper_TClassManip(TClass*);
   static void *new_FlowAnalysis_Helper(void *p = nullptr);
   static void *newArray_FlowAnalysis_Helper(Long_t size, void *p);
   static void delete_FlowAnalysis_Helper(void *p);
   static void deleteArray_FlowAnalysis_Helper(void *p);
   static void destruct_FlowAnalysis_Helper(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FlowAnalysis_Helper*)
   {
      ::FlowAnalysis_Helper *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::FlowAnalysis_Helper));
      static ::ROOT::TGenericClassInfo 
         instance("FlowAnalysis_Helper", "FlowAnalysis_Helper.h", 51,
                  typeid(::FlowAnalysis_Helper), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &FlowAnalysis_Helper_Dictionary, isa_proxy, 4,
                  sizeof(::FlowAnalysis_Helper) );
      instance.SetNew(&new_FlowAnalysis_Helper);
      instance.SetNewArray(&newArray_FlowAnalysis_Helper);
      instance.SetDelete(&delete_FlowAnalysis_Helper);
      instance.SetDeleteArray(&deleteArray_FlowAnalysis_Helper);
      instance.SetDestructor(&destruct_FlowAnalysis_Helper);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FlowAnalysis_Helper*)
   {
      return GenerateInitInstanceLocal(static_cast<::FlowAnalysis_Helper*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FlowAnalysis_Helper*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *FlowAnalysis_Helper_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::FlowAnalysis_Helper*>(nullptr))->GetClass();
      FlowAnalysis_Helper_TClassManip(theClass);
   return theClass;
   }

   static void FlowAnalysis_Helper_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/Users/cz263673/alice/Kit_Analysis/Flow/./FlowAnalysis_Helper.h");
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_FlowAnalysis_Helper(void *p) {
      return  p ? new(p) ::FlowAnalysis_Helper : new ::FlowAnalysis_Helper;
   }
   static void *newArray_FlowAnalysis_Helper(Long_t nElements, void *p) {
      return p ? new(p) ::FlowAnalysis_Helper[nElements] : new ::FlowAnalysis_Helper[nElements];
   }
   // Wrapper around operator delete
   static void delete_FlowAnalysis_Helper(void *p) {
      delete (static_cast<::FlowAnalysis_Helper*>(p));
   }
   static void deleteArray_FlowAnalysis_Helper(void *p) {
      delete [] (static_cast<::FlowAnalysis_Helper*>(p));
   }
   static void destruct_FlowAnalysis_Helper(void *p) {
      typedef ::FlowAnalysis_Helper current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FlowAnalysis_Helper

namespace {
  void TriggerDictionaryInitialization_FlowAnalysis_Helper_cxx_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./FlowAnalysis_Helper.cxx",
nullptr
    };
    static const char* includePaths[] = {
"/Users/cz263673/alice/sw/osx_arm64/ROOT/v6-32-06-alice1-local7/include",
"/Users/cz263673/alice/sw/osx_arm64/O2Physics/dev_muonRealign-local2/include",
"/Users/cz263673/alice/sw/osx_arm64/KFParticle/v1.1-5-local7/include",
"/Users/cz263673/alice/sw/osx_arm64/O2/dev_mchAlign-local1/include",
"/Users/cz263673/alice/sw/osx_arm64/O2/dev_mchAlign-local1/include/GPU",
"/Users/cz263673/alice/sw/osx_arm64/RapidJSON/v1.1.0-alice2-local1/include",
"/Users/cz263673/alice/sw/osx_arm64/ONNXRuntime/v1.18.1-local3/include/onnxruntime",
"/Users/cz263673/alice/sw/osx_arm64/FairMQ/v1.9.1-local1/include/fairmq",
"/Users/cz263673/alice/sw/osx_arm64/FairMQ/v1.9.1-local1/include",
"/Users/cz263673/alice/sw/osx_arm64/fastjet/v3.4.1_1.052-alice2-local3/include",
"/Users/cz263673/alice/sw/osx_arm64/JAliEn-ROOT/0.7.14-local8/include",
"/Users/cz263673/alice/sw/osx_arm64/ms_gsl/4.0.0-local1/include",
"/Users/cz263673/alice/sw/osx_arm64/Common-O2/v1.6.3-local3/include",
"/Users/cz263673/alice/sw/osx_arm64/Monitoring/v3.18.1-local3/include",
"/Users/cz263673/alice/sw/osx_arm64/libInfoLogger/v2.7.3-local3/include",
"/Users/cz263673/alice/sw/osx_arm64/HepMC3/3.3.0-local7/include",
"/Users/cz263673/alice/sw/osx_arm64/FairRoot/v18.4.9-alice3-local7/include",
"/Users/cz263673/alice/sw/osx_arm64/FairLogger/v1.11.1-local1/include",
"/Users/cz263673/alice/sw/osx_arm64/fmt/10.1.1-local1/include",
"/Users/cz263673/alice/sw/osx_arm64/GEANT3/v4-4-local7/include/TGeant3",
"/Users/cz263673/alice/sw/osx_arm64/GEANT4_VMC/v6-6-p3-local1/include/g4root",
"/Users/cz263673/alice/sw/osx_arm64/GEANT4_VMC/v6-6-p3-local1/include/geant4vmc",
"/Users/cz263673/alice/sw/osx_arm64/MCStepLogger/v0.6.1-local7/include",
"/Users/cz263673/alice/sw/osx_arm64/VMC/v2-0-local7/include/vmc",
"/Users/cz263673/alice/sw/osx_arm64/TBB/v2021.5.0-local1/include",
"/Users/cz263673/alice/sw/osx_arm64/Vc/1.4.1-local1/include",
"/Users/cz263673/alice/sw/osx_arm64/pythia/v8311-local3/include",
"/Users/cz263673/alice/sw/osx_arm64/boost/v1.83.0-alice2-local3/include",
"/Users/cz263673/alice/sw/osx_arm64/ROOT/v6-32-06-alice1-local7/etc/",
"/Users/cz263673/alice/sw/osx_arm64/ROOT/v6-32-06-alice1-local7/etc//cling",
"/Users/cz263673/alice/sw/osx_arm64/ROOT/v6-32-06-alice1-local7/etc//cling/plugins/include",
"/Users/cz263673/alice/sw/osx_arm64/ROOT/v6-32-06-alice1-local7/include/",
"/Users/cz263673/alice/sw/osx_arm64/ROOT/v6-32-06-alice1-local7/include",
"/Users/cz263673/alice/sw/osx_arm64/ROOT/v6-32-06-alice1-local7/include/",
"/Users/cz263673/alice/Kit_Analysis/Flow/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "FlowAnalysis_Helper_cxx_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$FlowAnalysis_Helper.h")))  __attribute__((annotate("$clingAutoload$./FlowAnalysis_Helper.cxx")))  FlowAnalysis_Helper;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "FlowAnalysis_Helper_cxx_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./FlowAnalysis_Helper.cxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"", payloadCode, "@",
"FlowAnalysis_Helper", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("FlowAnalysis_Helper_cxx_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_FlowAnalysis_Helper_cxx_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_FlowAnalysis_Helper_cxx_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_FlowAnalysis_Helper_cxx_ACLiC_dict() {
  TriggerDictionaryInitialization_FlowAnalysis_Helper_cxx_ACLiC_dict_Impl();
}
