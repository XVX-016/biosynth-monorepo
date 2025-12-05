
ALTER TABLE molecules 
ADD COLUMN IF NOT EXISTS molfile TEXT;

SELECT column_name, data_type, is_nullable 
FROM information_schema.columns 
WHERE table_name = 'molecules' AND column_name = 'molfile';


