// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME FlowAnalysis_Fitting_rflx
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

// Header files passed as explicit arguments
#include "FlowAnalysis_Fitting.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *FlowAnalysis_Fitting_Dictionary();
   static void FlowAnalysis_Fitting_TClassManip(TClass*);
   static void *new_FlowAnalysis_Fitting(void *p = nullptr);
   static void *newArray_FlowAnalysis_Fitting(Long_t size, void *p);
   static void delete_FlowAnalysis_Fitting(void *p);
   static void deleteArray_FlowAnalysis_Fitting(void *p);
   static void destruct_FlowAnalysis_Fitting(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FlowAnalysis_Fitting*)
   {
      ::FlowAnalysis_Fitting *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::FlowAnalysis_Fitting));
      static ::ROOT::TGenericClassInfo 
         instance("FlowAnalysis_Fitting", "", 34,
                  typeid(::FlowAnalysis_Fitting), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &FlowAnalysis_Fitting_Dictionary, isa_proxy, 0,
                  sizeof(::FlowAnalysis_Fitting) );
      instance.SetNew(&new_FlowAnalysis_Fitting);
      instance.SetNewArray(&newArray_FlowAnalysis_Fitting);
      instance.SetDelete(&delete_FlowAnalysis_Fitting);
      instance.SetDeleteArray(&deleteArray_FlowAnalysis_Fitting);
      instance.SetDestructor(&destruct_FlowAnalysis_Fitting);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FlowAnalysis_Fitting*)
   {
      return GenerateInitInstanceLocal(static_cast<::FlowAnalysis_Fitting*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::FlowAnalysis_Fitting*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *FlowAnalysis_Fitting_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::FlowAnalysis_Fitting*>(nullptr))->GetClass();
      FlowAnalysis_Fitting_TClassManip(theClass);
   return theClass;
   }

   static void FlowAnalysis_Fitting_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_FlowAnalysis_Fitting(void *p) {
      return  p ? new(p) ::FlowAnalysis_Fitting : new ::FlowAnalysis_Fitting;
   }
   static void *newArray_FlowAnalysis_Fitting(Long_t nElements, void *p) {
      return p ? new(p) ::FlowAnalysis_Fitting[nElements] : new ::FlowAnalysis_Fitting[nElements];
   }
   // Wrapper around operator delete
   static void delete_FlowAnalysis_Fitting(void *p) {
      delete (static_cast<::FlowAnalysis_Fitting*>(p));
   }
   static void deleteArray_FlowAnalysis_Fitting(void *p) {
      delete [] (static_cast<::FlowAnalysis_Fitting*>(p));
   }
   static void destruct_FlowAnalysis_Fitting(void *p) {
      typedef ::FlowAnalysis_Fitting current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::FlowAnalysis_Fitting

namespace {
  void TriggerDictionaryInitialization_FlowAnalysis_Fitting_rflx_Impl() {
    static const char* headers[] = {
"0",
nullptr
    };
    static const char* includePaths[] = {
"/Users/cz263673/alice/sw/osx_arm64/ROOT/v6-32-02-alice1-local3/include/",
"/Users/cz263673/alice/Kit_Analysis/Flow/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "FlowAnalysis_Fitting_rflx dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class FlowAnalysis_Fitting;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "FlowAnalysis_Fitting_rflx dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#ifndef FLOWANALYSIS_FITTING_H
#define FLOWANALYSIS_FITTING_H

#include <TCanvas.h>
#include <TH1F.h>
#include <TList.h>

#include <RooAbsData.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooChebychev.h>
#include <RooCrystalBall.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooHist.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooPowerSum.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

using namespace std;
using namespace RooFit;

enum ModelType { CB2 = 0, Chebychev, VWG, POL, Exp2, PolExp };

class FlowAnalysis_Fitting {
public:
  void init();
  void setModel(int flag_sig, int flag_bkg);
  void setChi2Max(double chi2) { mchi2max = chi2; };
  void runFitting(TH1D *hs, TList *ls, double ptmin, double ptmax,
                  double massmin, double massmax);
  void Print();

private:
  void CreateModel(RooWorkspace &w, RooRealVar x, int flag);
  int mflag_sig{0};
  int mflag_bkg{0};
  double mchi2max{1.};
};
#endif
#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"FlowAnalysis_Fitting", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("FlowAnalysis_Fitting_rflx",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_FlowAnalysis_Fitting_rflx_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_FlowAnalysis_Fitting_rflx_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_FlowAnalysis_Fitting_rflx() {
  TriggerDictionaryInitialization_FlowAnalysis_Fitting_rflx_Impl();
}
