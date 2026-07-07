#include "label.H"

namespace Foam
{
    class base
    {
        protected:
            label value;

        public:

            //- Construct from region mesh
            base(label value);

            //- Destructor
            virtual ~base()
            {}

            //- Disallow default bitwise copy construction
            base
            (
                const base&
            ) = delete;

            //- Disallow default bitwise assignment
            void operator=(const base&) = delete;

            //- test method
            virtual void method();

            //- second method
            virtual label second()
            {
                return 0.0;
            }
    };
}