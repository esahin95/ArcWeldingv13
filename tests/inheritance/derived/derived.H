#include "base.H"

namespace Foam
{
    class derived
    :
        public base
    {
        public:

            //- Construct from region mesh
            derived(label value);

            //- Destructor
            virtual ~derived()
            {}

            //- Disallow default bitwise copy construction
            derived
            (
                const derived&
            ) = delete;

            //- Disallow default bitwise assignment
            void operator=(const derived&) = delete;

            //- test method
            virtual void method();

            //- second method
            virtual label second()
            {
                return 1.0;
            }
    };
}