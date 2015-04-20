
#ifndef KRATOS_EXPORT_DLL_H
#define KRATOS_EXPORT_DLL_H

#ifdef KratosCore_BUILT_AS_STATIC
#  define KRATOS_EXPORT_DLL
#  define KRATOSCORE_NO_EXPORT
#else
#  ifndef KRATOS_EXPORT_DLL
#    ifdef KratosCore_EXPORTS
        /* We are building this library */
#      define KRATOS_EXPORT_DLL __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define KRATOS_EXPORT_DLL __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef KRATOSCORE_NO_EXPORT
#    define KRATOSCORE_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef KRATOSCORE_DEPRECATED
#  define KRATOSCORE_DEPRECATED __attribute__ ((__deprecated__))
#  define KRATOSCORE_DEPRECATED_EXPORT KRATOS_EXPORT_DLL __attribute__ ((__deprecated__))
#  define KRATOSCORE_DEPRECATED_NO_EXPORT KRATOSCORE_NO_EXPORT __attribute__ ((__deprecated__))
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define KRATOSCORE_NO_DEPRECATED
#endif

#endif
