#include "cache_base.h"
#include "utils.h"
#include "log.h"
#include "rng.h"
#include "address_home_lookup.h"

#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>

CacheBase::CacheBase(
   String name, UInt32 num_sets, UInt32 associativity, UInt32 cache_block_size,
   CacheBase::hash_t hash, AddressHomeLookup *ahl, int shared_cores, core_id_t core_id, bool is_last_level_cache, String partition_info)
:
   m_name(name),
   m_cache_size(UInt64(num_sets) * associativity * cache_block_size),
   m_associativity(associativity),
   m_blocksize(cache_block_size),
   m_hash(hash),
   m_num_sets(num_sets),
   m_ahl(ahl),
   m_shared_cores(shared_cores),
   m_core_id(core_id),
   m_is_last_level_cache(is_last_level_cache),
   m_partition_info(partition_info)

{
   m_log_blocksize = floorLog2(m_blocksize);
   m_log_num_sets = floorLog2(m_num_sets);

   is_hash_initialized = false;

   setStartArray = (int *) calloc(m_shared_cores + 1, sizeof(int));
   setLengthArray = (int *) calloc(m_shared_cores + 1, sizeof(int));

   // if(m_is_last_level_cache)
   //    std::cout << "LLC Partition Info: " << m_partition_info  << std::endl;

   LOG_ASSERT_ERROR((m_num_sets == (1UL << floorLog2(m_num_sets))) || (hash != CacheBase::HASH_MASK),
      "Caches of non-power of 2 size need funky hash function");
}

CacheBase::~CacheBase()
{}

// utilities
CacheBase::hash_t
CacheBase::parseAddressHash(String hash_name)
{
   if (hash_name == "mask")
      return CacheBase::HASH_MASK;
   else if (hash_name == "mod")
      return CacheBase::HASH_MOD;
   else if (hash_name == "rng1_mod")
      return CacheBase::HASH_RNG1_MOD;
   else if (hash_name == "rng2_mod")
      return CacheBase::HASH_RNG2_MOD;
   else if (hash_name == "prime_dis")
   	  return CacheBase::HASH_PRIME_DIS;
   else if (hash_name == "xor_mod")
   	  return CacheBase::HASH_XOR_MOD;
   else if (hash_name == "mersenne_mod")
   	  return CacheBase::HASH_MER_MOD;
   else if (hash_name == "fair_share")
   	  return CacheBase::HASH_FAIR_SHARE;
   else if (hash_name == "unfair_share")
   	  return CacheBase::HASH_UNFAIR_SHARE;
   else
      LOG_PRINT_ERROR("Invalid address hash function %s", hash_name.c_str());
}

void
CacheBase::splitAddress(const IntPtr addr, IntPtr& tag, UInt32& set_index) const
{
   tag = addr >> m_log_blocksize;

   IntPtr linearAddress = m_ahl ? m_ahl->getLinearAddress(addr) : addr;
   IntPtr block_num = linearAddress >> m_log_blocksize;

   switch(m_hash)
   {
      case CacheBase::HASH_MASK:
         set_index = block_num & (m_num_sets-1);
         break;
      case CacheBase::HASH_MOD:
         set_index = block_num % m_num_sets;
         break;
      case CacheBase::HASH_RNG1_MOD:
      {
         UInt64 state = rng_seed(block_num);
         set_index = rng_next(state) % m_num_sets;
         break;
      }
      case CacheBase::HASH_RNG2_MOD:
      {
         UInt64 state = rng_seed(block_num);
         rng_next(state);
         set_index = rng_next(state) % m_num_sets;
         break;
      }
      case CacheBase::HASH_PRIME_DIS:
      {

      	 //Prime Displacement Hashing Function (HPCA04)
         UInt64 si = block_num % m_num_sets;
         UInt64 Ti = block_num >> m_log_num_sets;
         UInt64 rho = 3;
         set_index = (rho * Ti + si) % m_num_sets;
         break;
      }
      case CacheBase::HASH_XOR_MOD:
      {
      	 //Based on related work of "Eliminating Conflict Misses Using Prime Number-Based Cache Indexing" (TC may 2005)
      	 //XOR based hash function
      	 UInt64 si = block_num % m_num_sets;
      	 UInt64 ti = (block_num >> m_log_num_sets) % m_num_sets;
      	 set_index = (si ^ ti); // ^ -> bitwise XOR
      	 break;
      }
      case CacheBase::HASH_MER_MOD:
      {
      	 //Based on related work of "Eliminating Conflict Misses Using Prime Number-Based Cache Indexing" (TC may 2005)
      	 //Mersenne based hash function
      	 //Disadvantage of this mod is that we will not use one set in our cache.
      	 set_index = block_num % (m_num_sets - 1);
      	 break;
      }

      default:
         LOG_PRINT_ERROR("Invalid hash function %d", m_hash);
   }
}

// Only the LLC supports Fair and unfair partitioning scheme
void
CacheBase::splitAddress(const core_id_t core_id, const IntPtr addr, IntPtr& tag, UInt32& set_index) const
{
   tag = addr >> m_log_blocksize;

   IntPtr linearAddress = m_ahl ? m_ahl->getLinearAddress(addr) : addr;
   IntPtr block_num = linearAddress >> m_log_blocksize;

   switch(m_hash)
   {
      case CacheBase::HASH_MASK:
         set_index = block_num & (m_num_sets-1);
         break;
      case CacheBase::HASH_MOD:
         set_index = block_num % m_num_sets;
         break;
      case CacheBase::HASH_RNG1_MOD:
      {
         UInt64 state = rng_seed(block_num);
         set_index = rng_next(state) % m_num_sets;
         break;
      }
      case CacheBase::HASH_RNG2_MOD:
      {
         UInt64 state = rng_seed(block_num);
         rng_next(state);
         set_index = rng_next(state) % m_num_sets;
         break;
      }
      case CacheBase::HASH_PRIME_DIS:
      {

      	 //Prime Displacement Hashing Function (HPCA04)
         UInt64 si = block_num % m_num_sets;
         UInt64 Ti = block_num >> m_log_num_sets;
         UInt64 rho = 3;
         set_index = (rho * Ti + si) % m_num_sets;
         break;
      }
      case CacheBase::HASH_XOR_MOD:
      {
      	 //Based on related work of "Eliminating Conflict Misses Using Prime Number-Based Cache Indexing" (TC may 2005)
      	 //XOR based hash function
      	 UInt64 si = block_num % m_num_sets;
      	 UInt64 ti = (block_num >> m_log_num_sets) % m_num_sets;
      	 set_index = (si ^ ti); // ^ -> bitwise XOR
      	 break;
      }
      case CacheBase::HASH_MER_MOD:
      {
      	 //Based on related work of "Eliminating Conflict Misses Using Prime Number-Based Cache Indexing" (TC may 2005)
      	 //Mersenne based hash function
      	 //Disadvantage of this mod is that we will not use one set in our cache.
      	 set_index = block_num % (m_num_sets - 1);
      	 break;
      }
      case CacheBase::HASH_FAIR_SHARE:
      {
         /* For LLC, in Set-Based Partitioning
            Each core is given equal number of
            Sets in the cache.
            For non-LLC, it is same as HASH_MOD
          */
         set_index = ((block_num << floorLog2(m_shared_cores)) + core_id) % (m_num_sets);
         break;
      }

      case CacheBase::HASH_UNFAIR_SHARE:
      {
         initializeUnfairHash();

         if(isAddressShared(addr)){
            LOG_PRINT("Access Shared Address: 0x%x, core: %d", addr, core_id);
            set_index = setStartArray[4] + block_num % setLengthArray[4];
         }
         else{
            set_index = setStartArray[core_id] + block_num % setLengthArray[core_id];
            LOG_PRINT("Access Private Address: 0x%x, core: %d", addr, core_id);
         }

         break;
      }

      default:
         LOG_PRINT_ERROR("Invalid hash function %d", m_hash);
   }
}

int* parsePartitionInfo(String partitionInfoS){
   int numPartitions = 1;

   int *coreShare = NULL;

   std::string partitionInfo = std::string(partitionInfoS.c_str());

   std::string::iterator it;
   for(it = partitionInfo.begin(); it != partitionInfo.end(); it++){
      if(*it == ':')
         numPartitions++;
   }

   coreShare = (int *) malloc(sizeof(int) * numPartitions);

   std::stringstream ss(partitionInfo);

   std::string s;
   int i = 0;
   while (std::getline(ss, s, ':')) {
      coreShare[i] = std::stoi(s);
      i++;
   }

   return coreShare;
}

bool CacheBase::isAddressShared(const IntPtr addr) const{
   printf("Address Access: 0x%lx\n", addr);
   IntPtr temp_addr;
   for(int i = 0; i < 79; i++){
      temp_addr = 0xfffffffff & addr;
      if(temp_addr >= sharedAddressRanges[i][0] && temp_addr <= sharedAddressRanges[i][1])
         return true;
   }

   return false;
}

void CacheBase::initializeUnfairHash() const{
   if(is_hash_initialized)
      return;

   FILE *filePointer;
   
   if ((filePointer = fopen("/media/madhur/CommonSpace/Work/SystemSimulators/Sniper/snipersim/test/tracing/common_ranges.txt","r")) == NULL){
       printf("Error! opening file");
       // Program exits if the file pointer returns NULL.
       exit(1);
   }

   for(int i = 0; i < 79; i++){
      fscanf(filePointer,"{%lx, %lx}\n", &sharedAddressRanges[i][0], &sharedAddressRanges[i][1]);
   }

   // for(int i = 0; i < 962; i++){
   //    printf("{%lx, %lx}\n", sharedAddressRanges[i][0], sharedAddressRanges[i][1]);
   // }

   is_hash_initialized = true;

   int *coreShare = NULL;

   coreShare = parsePartitionInfo(m_partition_info);
   
   sharedAddressStart   =  0x000000000000;
   sharedAddressEnd     =  0xffffffffffff;

   std::cout << "Number Of Sets: " << m_num_sets << std::endl;

   for(int i = 0; i < m_shared_cores + 1; i++){
      setLengthArray[i] = (int)((coreShare[i] / 100.0) * m_num_sets);

      if(i == 0){
         setStartArray[i] = 0;
      } else{
         setStartArray[i] = setStartArray[i - 1] + setLengthArray[i - 1];
      }

      std::cout << "Core " << i << ":\tSet Start: " << setStartArray[i] << "\tSet Length: " << setLengthArray[i] << std::endl;
   }
}

void
CacheBase::splitAddress(const IntPtr addr, IntPtr& tag, UInt32& set_index,
                  UInt32& block_offset) const
{
   block_offset = addr & (m_blocksize-1);
   splitAddress(addr, tag, set_index);
}

void
CacheBase::splitAddress(const core_id_t core_id, const IntPtr addr, IntPtr& tag, UInt32& set_index,
                  UInt32& block_offset) const
{
   block_offset = addr & (m_blocksize-1);
   splitAddress(core_id, addr, tag, set_index);
}

IntPtr
CacheBase::tagToAddress(const IntPtr tag)
{
   return tag << m_log_blocksize;
}
